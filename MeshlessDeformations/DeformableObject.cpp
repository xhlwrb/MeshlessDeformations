#include "DeformableObject.h"
#include "Globals.h"
#include "matrix/newmat.h"
#include "matrix/newmatap.h"
#include "matrix/tmt.h"

using namespace Model;

//CONSTRUCTORS
DeformableObject::DeformableObject(core::stringw name) : ModelObject(ObjectType::DeformableObject)
{
	//beta 
	// To control the amount of influence the linear transformation has on the shape matching, the parameter beta is used
	// By using a smaller value for β, the ratio of rotation-to-pure deformation will increase
	// and so the goal positions will more closely match the original, undeformed shape of the object.
	beta = 0.0f;
	//alpha 
	// a parameter which simulates stiffness
	alpha = 0.5f;
	elasticity = 0.5f;
	this->name = name;
	finished = false;
	isSelected = false;
}

void DeformableObject::doInitialCalculations() {
	// calculate centre of mass (t_0)
	vector3df tmpSum(0.0f,0.0f,0.0f);
	float tmpMass = 0.0f;
	for (u32 index = 0; index < particles.size(); index++) {
		Particle* particle = particles[index];
		//tmpSum += particle->originalPosition * particle->weight;
		// sigma(m_i * x0_i)
		tmpSum += particle->weight * particle->originalPosition;
		tmpMass += particle->weight;
	}
	// x0_cm = sigma(m_i * x0_i) / sigma(m_i)
	originalCentreOfMass = tmpSum / tmpMass;   

	// update q_i
	//clear()
	// Clears the array and deletes all allocated memory.
	q.clear();
	for (u32 index = 0; index < particles.size(); index++) {
		Particle* p = particles[index];
		//push_back()
		// Adds an element at back of array.
		// If the array is too small to add this new element it is made bigger.
		// q_i = x0_i - x0_cm
		q.push_back(p->originalPosition - originalCentreOfMass);
	}

	// calculate A_qq
	matrix4 A_qq_inverse;
	for(int i = 0; i < (int)particles.size(); i++) {
		// m_i * q_i * (q_i)^T
		Particle *par = particles[i];
		f32 m = par->weight;
		A_qq_inverse(0,0)+=m*q[i].X*q[i].X; A_qq_inverse(1,0)+=m*q[i].X*q[i].Y; A_qq_inverse(2,0)+=m*q[i].X*q[i].Z;
		A_qq_inverse(0,1)+=m*q[i].Y*q[i].X; A_qq_inverse(1,1)+=m*q[i].Y*q[i].Y; A_qq_inverse(2,1)+=m*q[i].Y*q[i].Z;
		A_qq_inverse(0,2)+=m*q[i].Z*q[i].X; A_qq_inverse(1,2)+=m*q[i].Z*q[i].Y; A_qq_inverse(2,2)+=m*q[i].Z*q[i].Z;
	}
	// ^T
	A_qq_inverse.getInverse(A_qq);
	// update q_tilde_i
	// q_tilde = [q_x, q_y, q_z, q_x^2, q_y^2, q_z^2, q_x * q_y, q_y * q_z, q_z * q_x]^T
	q_tilde.clear();
	for (u32 index = 0; index < particles.size(); index++) {
		Particle* p = particles[index];
		ColumnVector t(9);
		t.element(0) = q[index].X;
		t.element(1) = q[index].Y;
		t.element(2) = q[index].Z;
		t.element(3) = q[index].X * q[index].X;
		t.element(4) = q[index].Y * q[index].Y;
		t.element(5) = q[index].Z * q[index].Z;
		t.element(6) = q[index].X * q[index].Y;
		t.element(7) = q[index].Y * q[index].Z;
		t.element(8) = q[index].Z * q[index].X;
		q_tilde.push_back(t);
	}

	// calculate A_tilde_qq
	SquareMatrix A_tilde_qq_inverse(9);
	A_tilde_qq_inverse = 0.0;
	for (u32 index = 0; index < particles.size(); index++) {
		// m_i * q_tilde_i * q_tilde_i^T
		A_tilde_qq_inverse += q_tilde[index] * q_tilde[index].t();
	}

	A_tilde_qq = SquareMatrix(9);
	A_tilde_qq = 0.0;
	try {
        A_tilde_qq = A_tilde_qq_inverse.i();
	}catch (SingularException) {
		A_tilde_qq = IdentityMatrix(9);
	}
}

void DeformableObject::addParticle(Particle* particle) {
	int nr = particles.size();
	particles.push_back( particle );
	particle->nr = nr;
	finished = false;
}

array<Particle*> DeformableObject::getParticles(void) {
	return particles;
}

matrix4 DeformableObject::calculateA_pqMatrix()
{
	matrix4 A_pq;
	for(int i=0;i<(int)particles.size(); i++) {
		// m_i * p_i * (q_i)^T
		Particle *par = particles[i];
		f32 m = par->weight;
		A_pq(0,0)+=m*p[i].X*q[i].X; A_pq(1,0)+=m*p[i].X*q[i].Y; A_pq(2,0)+=m*p[i].X*q[i].Z;
		A_pq(0,1)+=m*p[i].Y*q[i].X; A_pq(1,1)+=m*p[i].Y*q[i].Y; A_pq(2,1)+=m*p[i].Y*q[i].Z;
		A_pq(0,2)+=m*p[i].Z*q[i].X; A_pq(1,2)+=m*p[i].Z*q[i].Y; A_pq(2,2)+=m*p[i].Z*q[i].Z;
	}
//██████████████████████████████████████████████████████████████████████████████████
	return A_pq.getTransposed();
}

matrix4 DeformableObject::calculateRotationMatrix(matrix4 A_pq)
{
	try {
		// sqrt() wrapped for matrix4 in globals.cpp
		matrix4 S = sqrt(A_pq.getTransposed() * A_pq);

		matrix4 Sinverse;
		if(S.getInverse(Sinverse) == false)
			std::cout << "no inverse for S" << std::endl;

		return A_pq * Sinverse;

	} catch(ConvergenceException) {
		return matrix4();
	}
}

void DeformableObject::updateGoalPositions(f32 timeElapsed) {

	// calc roofPosition only based on external forces 
	for (u32 index = 0; index < particles.size(); index++) {
		Particle* par = particles[index];
		par->velocityRoof = par->velocity + par->getForces() * timeElapsed;
		par->roofPosition = par->position + par->velocityRoof * timeElapsed ;
	}	

	vector3df tmpSum(0.0f,0.0f,0.0f);
	float tmpMass = 0.0f;
	for (u32 index = 0; index < particles.size(); index++) {
		Particle* par = particles[index];
		tmpSum += par->roofPosition * par->weight;
		tmpMass += par->weight;
	}
	currentCentreOfMass = tmpSum / tmpMass;

	// update p_i
	p.clear();
	for (u32 index = 0; index < particles.size(); index++) {
		Particle* par = particles[index];
		p.push_back(par->roofPosition - currentCentreOfMass);
	}

	if (Globals::mode == DeformationMode::Basic) {
		// 1. Rigid Body Dynamics 
		matrix4 A_pq = calculateA_pqMatrix();
		matrix4 R = calculateRotationMatrix(A_pq);

		// calculate goal positions and integrate
		for (u32 i = 0; i < particles.size(); i++) {
			vector3df goalPosition = q[i];
			R.rotateVect(goalPosition);
			goalPosition = goalPosition + currentCentreOfMass;

			if(particles[i]->canMoveTo(goalPosition)) { 
				particles[i]->goalPosition = goalPosition;
			}
		}
	} else if (Globals::mode == DeformationMode::Linear) {
		// 2. Linear Deformations
		matrix4 A_pq = calculateA_pqMatrix();
		matrix4 R = calculateRotationMatrix(A_pq);
//██████████████████████████████████████████████████████████████████████████████████
		matrix4 A = volumeNormalize(A_pq * A_qq);
		matrix4 Transform;
		for (int i=0; i<16; i++) {
			A.M[i] *= beta;
			R.M[i] *= (1-beta);
			Transform.M[i] = R.M[i] + A.M[i];
		}

		// calculate goal positions and integrate
		for (u32 i = 0; i < particles.size(); i++) {
			vector3df goalPosition = q[i];
			Transform.rotateVect(goalPosition);
			goalPosition = goalPosition + currentCentreOfMass;

			if(particles[i]->canMoveTo(goalPosition)) { 
				particles[i]->goalPosition = goalPosition;
			}
		}
	} else if (Globals::mode == DeformationMode::Quadratic) {
		// 3. Quadratic deformation
		Matrix T_tilde = Matrix(3,9); 
		T_tilde = 0.0;

		p_tilde.clear();
		for (u32 index = 0; index < particles.size(); index++) {
			Particle* par = particles[index];
			ColumnVector t(3);
			t.element(0) = par->roofPosition.X - currentCentreOfMass.X;
			t.element(1) = par->roofPosition.Y - currentCentreOfMass.Y;
			t.element(2) = par->roofPosition.Z - currentCentreOfMass.Z;
			p_tilde.push_back(t);
		}
		Matrix A_tilde_pq(3, 9);
		A_tilde_pq = 0.0;
		for (u32 index = 0; index < particles.size(); index++) {
			A_tilde_pq += p_tilde[index] * q_tilde[index].t();
		}
		Matrix A_tilde(3, 9);
		A_tilde = 0.0;
		A_tilde = (A_tilde_pq * A_tilde_qq);


		Matrix A_tilde_square(9, 9);
//██████████████████████████████████████████████████████████████████████████████████
		A_tilde_square = IdentityMatrix(9);
		for(int j=0;j<3;j++) {
			for(int i=0;i<9;i++) { 
				A_tilde_square.element(j, i) = A_tilde.element(j, i);
			}
		}
		double det = Determinant(A_tilde_square);
		double cbrt = pow( fabs(det), 1.0/9.0 );
		cbrt = ( det < 0 ) ? -cbrt : cbrt;
		A_tilde = A_tilde / cbrt;

		matrix4 A_pq = calculateA_pqMatrix();
		matrix4 R = calculateRotationMatrix(A_pq);
		Matrix R_tilde = Matrix(3,9);
		R_tilde = 0.0;
		for(int j=0;j<3;j++) {
			for(int i=0;i<3;i++) { 
				R_tilde.element(j, i) = R(j, i);
			}
		}
		matrix4 A = volumeNormalize(A_pq * A_qq);
		T_tilde = beta * A_tilde + (1-beta) * R_tilde;

		// calculate goal positions and integrate
		for (u32 i = 0; i < particles.size(); i++) {
			ColumnVector goalN(3);
			goalN = T_tilde * q_tilde[i];
			vector3df goalPosition = vector3df(goalN.element(0), goalN.element(1), goalN.element(2))  + currentCentreOfMass;
			
			if(particles[i]->canMoveTo(goalPosition)) { 
				particles[i]->goalPosition = goalPosition;
			}
		}
	}
}

void DeformableObject::update(f32 timeElapsed) {

	if(!finished) {
		doInitialCalculations();
		finished = true;
	}

	updateGoalPositions(timeElapsed);

	for (u32 i = 0; i < particles.size(); i++) {
		Particle* particle = particles [i];
		if (!Globals::stopMode) {
			particle->update(timeElapsed);
		}
	}
}

void DeformableObject::setVisible(bool isVisible){
	//setVisible
	// Sets if the node should be visible or not.
	ISceneNode::setVisible(isVisible);

	for (u32 i = 0; i < particles.size(); i++) {
		Particle* particle = particles [i];
		if (isVisible) {
//██████████████████████████████████████████████████████████████████████████████████
			particle->setID(particle->getID() | ObjectType::SelectableObject);
		} else {
			particle->setID(particle->getID() ^ ObjectType::SelectableObject);
		}
	}
		
}

/*
#Code shared across examples
import pylab, string, sys

def stdDev(X):
    mean = sum(X)/float(len(X))
    tot = 0.0
    for x in X:
        tot += (x - mean)**2
    return (tot/len(X))**0.5

def scaleFeatures(vals):
    """Assumes vals is a sequence of numbers"""
    pylab.plot(vals)
    result = pylab.array(vals)
    mean = sum(result)/float(len(result))
    result = result - mean
    sd = stdDev(result)
    result = result/sd
    pylab.plot(result)
    pylab.show()
    return result

class Point(object):
    def __init__(self, name, originalAttrs):
        """originalAttrs is an array"""
        self.name = name
        self.attrs = originalAttrs
    def dimensionality(self):
        return len(self.attrs)
    def getAttrs(self):
        return self.attrs
    def distance(self, other):
        #Euclidean distance metric
        result = 0.0
        for i in range(self.dimensionality()):
            result += (self.attrs[i] - other.attrs[i])**2
        return result**0.5
    def getName(self):
        return self.name
    def toStr(self):
        return self.name + str(self.attrs)
    def __str__(self):
        return self.name        
    
class Cluster(object):
    """ A Cluster is defines as a set of elements, all having 
    a particular type """
    def __init__(self, points, pointType):
        """ Elements of a cluster are saved in self.points
        and the pointType is also saved """
        self.points = points
        self.pointType = pointType
    def singleLinkageDist(self, other):
        """ Returns the float distance between the points that 
        are closest to each other, where one point is from 
        self and the other point is from other. Uses the 
        Euclidean dist between 2 points, defined in Point."""
        dis = self.points[0].distance(other.points[0])
        for i in self.points:
            for j in other.points:
                disN = i.distance(j)
                if disN < dis:
                    dis = disN
        return dis
    def maxLinkageDist(self, other):
        """ Returns the float distance between the points that 
        are farthest from each other, where one point is from 
        self and the other point is from other. Uses the 
        Euclidean dist between 2 points, defined in Point."""
        dis = 0
        for i in self.points:
            for j in other.points:
                disN = i.distance(j)
                if disN > dis:
                    dis = disN
        return dis
    def averageLinkageDist(self, other):
        """ Returns the float average (mean) distance between all 
        pairs of points, where one point is from self and the 
        other point is from other. Uses the Euclidean dist 
        between 2 points, defined in Point."""
        dis = 0
        index = 0
        for i in self.points:
            for j in other.points:
                dis += i.distance(j)
                index += 1
        return float(dis) / index
    def members(self):
        for p in self.points:
            yield p
    def isIn(self, name):
        """ Returns True is the element named name is in the cluster
        and False otherwise """
        for p in self.points:
            if p.getName() == name:
                return True
        return False
    def toStr(self):
        result = ''
        for p in self.points:
            result = result + p.toStr() + ', '
        return result[:-2]
    def getNames(self):
        """ For consistency, returns a sorted list of all 
        elements in the cluster """
        names = []
        for p in self.points:
            names.append(p.getName())
        return sorted(names)
    def __str__(self):
        names = self.getNames()
        result = ''
        for p in names:
            result = result + p + ', '
        return result[:-2]

class ClusterSet(object):
    """ A ClusterSet is defined as a list of clusters """
    def __init__(self, pointType):
        """ Initialize an empty set, without any clusters """
        self.members = []
        self.pointType = pointType
    def add(self, c):
        """ Append a cluster to the end of the cluster list
        only if it doesn't already exist. If it is already in the 
        cluster set, raise a ValueError """
        if c in self.members:
            raise ValueError
        self.members.append(c)
    def getClusters(self):
        return self.members[:]
    def mergeClusters(self, c1, c2):
        """ Assumes clusters c1 and c2 are in self
        Adds to self a cluster containing the union of c1 and c2
        and removes c1 and c2 from self """
        # TO DO
        points = []
        type = self.pointType
        for point1 in c1.members():
            points.append(point1)
        for point2 in c2.members():
            if not c1.isIn(point2.getName()):
                points.append(point2)
        cN = Cluster(points, type)
        self.add(cN)
        self.members.remove(c1)
        self.members.remove(c2)
    def findClosest(self, linkage):
        """ Returns a tuple containing the two most similar 
        clusters in self
        Closest defined using the metric linkage """
        # TO DO
        min = None
        cluster1 = None
        cluster2 = None
        for i in range(len(self.getClusters())):
            for j in range(i+1, len(self.getClusters())):
                dis = linkage(self.getClusters()[i], self.getClusters()[j])
                if min == None or min > dis:
                    min = dis
                    cluster1 = self.getClusters()[i]
                    cluster2 = self.getClusters()[j]
        return cluster1, cluster2
    def mergeOne(self, linkage):
        """ Merges the two most simililar clusters in self
        Similar defined using the metric linkage
        Returns the clusters that were merged """
        # TO DO
        cluster1, cluster2 = self.findClosest(linkage)
        self.mergeClusters(cluster1, cluster2)
        return cluster1, cluster2
    def numClusters(self):
        return len(self.members)
    def toStr(self):
        cNames = []
        for c in self.members:
            cNames.append(c.getNames())
        cNames.sort()
        result = ''
        for i in range(len(cNames)):
            names = ''
            for n in cNames[i]:
                names += n + ', '
            names = names[:-2]
            result += '  C' + str(i) + ':' + names + '\n'
        return result

#City climate example
class City(Point):
    pass

def readCityData(fName, scale = False):
    """Assumes scale is a Boolean.  If True, features are scaled"""
    dataFile = open(fName, 'r')
    numFeatures = 0
    #Process lines at top of file
    for line in dataFile: #Find number of features
        if line[0:4] == '#end': #indicates end of features
            break
        numFeatures += 1
    numFeatures -= 1
    featureVals = []
    
    #Produce featureVals, cityNames
    featureVals, cityNames = [], []
    for i in range(numFeatures):
        featureVals.append([])
        
    #Continue processing lines in file, starting after comments
    for line in dataFile:
        dataLine = string.split(line[:-1], ',') #remove newline; then split
        cityNames.append(dataLine[0])
        for i in range(numFeatures):
            featureVals[i].append(float(dataLine[i+1]))
            
    #Use featureVals to build list containing the feature vectors
    #For each city scale features, if needed
    if scale:
        for i in range(numFeatures):
            featureVals[i] = scaleFeatures(featureVals[i])
    featureVectorList = []
    for city in range(len(cityNames)):
        featureVector = []
        for feature in range(numFeatures):
            featureVector.append(featureVals[feature][city])
        featureVectorList.append(featureVector)
    return cityNames, featureVectorList

def buildCityPoints(fName, scaling):
    cityNames, featureList = readCityData(fName, scaling)
    points = []
    for i in range(len(cityNames)):
        point = City(cityNames[i], pylab.array(featureList[i]))
        points.append(point)
    return points

#Use hierarchical clustering for cities
def hCluster(points, linkage, numClusters, printHistory):
    cS = ClusterSet(City)
    for p in points:
        cS.add(Cluster([p], City))
    history = []
    while cS.numClusters() > numClusters:
        merged = cS.mergeOne(linkage)
        history.append(merged)
    if printHistory:
        print ''
        for i in range(len(history)):
            names1 = []
            for p in history[i][0].members():
                names1.append(p.getName())
            names2 = []
            for p in history[i][1].members():
                names2.append(p.getName())
            print 'Step', i, 'Merged', names1, 'with', names2
            print ''
    print 'Final set of clusters:'
    print cS.toStr()
    return cS

def test():
    #points = buildCityPoints('/Users/Rose/Documents/University/Course Online/Introduction to Computational Thinking and Data Science/ProSet6/cityTemps.txt', False)
    #hCluster(points, Cluster.singleLinkageDist, 10, False)
    #hCluster(points, Cluster.singleLinkageDist, 5, False)
    points = buildCityPoints('/Users/Rose/Documents/University/Course Online/Introduction to Computational Thinking and Data Science/ProSet6/cityTemps.txt', True)
    hCluster(points, Cluster.singleLinkageDist, 5, False)
    #hCluster(points, Cluster.averageLinkageDist, 10, False)
    #hCluster(points, Cluster.singleLinkageDist, 10, False)

test()

*/
