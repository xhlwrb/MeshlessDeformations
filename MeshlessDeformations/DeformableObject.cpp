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

double DeformableObject::averageLinkageDist(core::array<Particle*> cluster1, core::array<Particle*> cluster2){
	/* Returns the float average (mean) distance between all
	pairs of points, where one point is from self and the
	other point is from other.Uses the Euclidean dist
	between 2 points, defined in Point.*/
	double dist = 0;
	int index = 0;
	for (int i = 0; i < cluster1.size(); i++) {
		for (int j = 0; j < cluster2.size(); j++){
			double diff_x = cluster1[i]->getPosition().X - cluster1[i]->getPosition().X;
			double diff_y = cluster1[i]->getPosition().Y - cluster1[i]->getPosition().Y;
			double diff_z = cluster1[i]->getPosition().Z - cluster1[i]->getPosition().Z;
			dist = sqrt(diff_x * diff_x + diff_y * diff_y + diff_z * diff_z);
			index += 1;
		}
	}
	double ave_dist = dist / index;
	return ave_dist;
}

core::array<core::array<Particle*>> DeformableObject::findClosest(core::array<core::array<Particle*>> clusters){
	double min = std::numeric_limits<double>::max();
	double dist = 0;
	core::array<Particle*> cluster1;
	core::array<Particle*> cluster2;
	for (int i = 0; i < clusters.size(); i++) {
		for (int j = i + 1; j < clusters.size(); j++) {
			dist = DeformableObject::averageLinkageDist(clusters[i], clusters[j]);
			if (min > dist){
				min = dist;
				cluster1 = clusters[i];
				cluster2 = clusters[j];
			}
		}
	}
	core::array<core::array<Particle*>> ret;
	ret.push_back(cluster1);
	ret.push_back(cluster2);
	return ret;
}

int DeformableObject::search_particle(core::array<Particle*>cluster, Particle* particle){
	for (int i = 0; i < cluster.size(); i++){
		if (cluster[i]->getPosition() == particle->getPosition()){
			return i;
		}
	}
	return -1;
}

int DeformableObject::search_cluster(core::array<core::array<Particle*>> clusters, core::array<Particle*>cluster){
	for (int i = 0; i < clusters.size(); i++){
		for (int j = 0; j < clusters[i].size(); j++){
			if (search_particle(cluster, clusters[i][j]) == -1){
				return -1;
			}
		}
		return i;
	}
}

core::array<core::array <Particle*>> DeformableObject::mergeClusters(core::array<core::array<Particle*>> clusters, core::array<Particle*> cluster1, core::array<Particle*> cluster2){
	core::array<Particle*> new_cluster;
	int inx_cluster1, inx_cluster2;
	for (int i = 0; i < cluster1.size(); i++) {
		new_cluster.push_back(cluster1[i]);
	}
	for (int i = 0; i < cluster2.size(); i++) {
			if (search_particle(cluster1, cluster2[i]) == -1){
			new_cluster.push_back(cluster2[i]);
		}
	}
	clusters.push_back(new_cluster);
	inx_cluster1 = search_cluster(clusters, cluster1);
	inx_cluster2 = search_cluster(clusters, cluster2);
	clusters.erase(inx_cluster1);
	clusters.erase(inx_cluster2);
	return clusters;
}

core::array<core::array<Particle*>> DeformableObject::mergeOne(core::array<core::array<Particle*>> clusters){
	core::array<Particle*> cluster1;
	core::array<Particle*> cluster2;
	core::array<core::array<Particle*>> ret_clusters;
	core::array<core::array<Particle*>> ret_findclosest = DeformableObject::findClosest(clusters);
	cluster1 = ret_findclosest[0];
	cluster2 = ret_findclosest[1];
	ret_clusters = DeformableObject::mergeClusters(clusters, cluster1, cluster2);
	return ret_clusters;
}

core::array<core::array<Particle*>> DeformableObject::cluster_merge(array<Particle*> all_particles, int num_clusters){
	core::array<core::array <Particle*>> clusters;
	core::array <Particle*> push;
	for (int i = 0; i < all_particles.size(); i++){
		push.push_back(all_particles[i]);
		clusters.push_back(push);
		push.clear();
	}
	while(clusters.size() > num_clusters){
		clusters = DeformableObject::mergeOne(clusters);
	}
	return clusters;
}

int DeformableObject::which_cluster(Particle* p){
	for (int i = 0; i < clusters.size(); i++){
		// this cluster
		if (search_particle(clusters[i], p) != -1){
			// in this cluster
			return i;
		}
	}
	return -1;
}

void DeformableObject::doInitialCalculations() {
	if (Globals::mode == DeformationMode::ClusterBased) {
		// 4. Cluster-Bases deformation

		// merge clusters
		int num_clusters = 5;
		DeformableObject::clusters = DeformableObject::cluster_merge(particles, num_clusters);

		// calculate centre of mass (t_0)
		core::array<vector3df> originalCentreOfMass, tmpSum;
		core::array<float> tmpMass = 0.0f;
		tmpSum = (0.0f, 0.0f, 0.0f);
		for (int i = 0; i < clusters.size(); i++){
			// this cluster
			vector3df sum(0.0f, 0.0f, 0.0f);
			float mass = 0.0f;
			for (u32 index = 0; index < particles.size(); index++) {
				Particle* particle = particles[index];
				//tmpSum += particle->originalPosition * particle->weight;
				// sigma(m_i * x0_i)
				if (which_cluster(particle) == i){
					// if this particle is in this cluster
					sum.X += particle->weight * particle->originalPosition.X;
					sum.Y += particle->weight * particle->originalPosition.Y;
					sum.Z += particle->weight * particle->originalPosition.Z;
					mass += particle->weight;
				}// this particle over
			}// all the particles in this cluster over
			tmpSum.push_back(sum);
			tmpMass.push_back(mass);
			// next cluster
		}// all the clusters over
		// x0_cm = sigma(m_i * x0_i) / sigma(m_i)
		for (int i = 0; i < clusters.size(); i++){
			vector3df centre(0.0f, 0.0f, 0.0f);
			centre.X += tmpSum[i].X / tmpMass[i];
			centre.Y += tmpSum[i].Y / tmpMass[i];
			centre.Z += tmpSum[i].Z / tmpMass[i];
			originalCentreOfMass.push_back(centre);
		}

		// update q_i
		//clear()
		// Clears the array and deletes all allocated memory.
		q.clear();
		for (u32 index = 0; index < particles.size(); index++) {
			Particle* particle = particles[index];
			//push_back()
			// Adds an element at back of array.
			// If the array is too small to add this new element it is made bigger.
			// q_i = x0_i - x0_cm
			int i = which_cluster(particle);
			if (i != -1){
				q.push_back(particle->originalPosition - originalCentreOfMass[i]);
			}
		}

		// calculate A_qq
		matrix4 A_qq_inverse;
		for (int i = 0; i < (int)particles.size(); i++) {
			// m_i * q_i * (q_i)^T
			Particle *par = particles[i];
			f32 m = par->weight;
			A_qq_inverse(0, 0) += m*q[i].X*q[i].X; A_qq_inverse(1, 0) += m*q[i].X*q[i].Y; A_qq_inverse(2, 0) += m*q[i].X*q[i].Z;
			A_qq_inverse(0, 1) += m*q[i].Y*q[i].X; A_qq_inverse(1, 1) += m*q[i].Y*q[i].Y; A_qq_inverse(2, 1) += m*q[i].Y*q[i].Z;
			A_qq_inverse(0, 2) += m*q[i].Z*q[i].X; A_qq_inverse(1, 2) += m*q[i].Z*q[i].Y; A_qq_inverse(2, 2) += m*q[i].Z*q[i].Z;
		}
		// ^T
		A_qq_inverse.getInverse(A_qq);
	}
	else{
		// calculate centre of mass (t_0)
		vector3df tmpSum(0.0f, 0.0f, 0.0f);
		float tmpMass = 0.0f;
		for (u32 index = 0; index < particles.size(); index++) {
			Particle* particle = particles[index];
			//tmpSum += particle->originalPosition * particle->weight;
			// sigma(m_i * x0_i)
			tmpSum.X += particle->weight * particle->originalPosition.X;
			tmpSum.Y += particle->weight * particle->originalPosition.Y;
			tmpSum.Z += particle->weight * particle->originalPosition.Z;
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
		for (int i = 0; i < (int)particles.size(); i++) {
			// m_i * q_i * (q_i)^T
			Particle *par = particles[i];
			f32 m = par->weight;
			A_qq_inverse(0, 0) += m*q[i].X*q[i].X; A_qq_inverse(1, 0) += m*q[i].X*q[i].Y; A_qq_inverse(2, 0) += m*q[i].X*q[i].Z;
			A_qq_inverse(0, 1) += m*q[i].Y*q[i].X; A_qq_inverse(1, 1) += m*q[i].Y*q[i].Y; A_qq_inverse(2, 1) += m*q[i].Y*q[i].Z;
			A_qq_inverse(0, 2) += m*q[i].Z*q[i].X; A_qq_inverse(1, 2) += m*q[i].Z*q[i].Y; A_qq_inverse(2, 2) += m*q[i].Z*q[i].Z;
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
		}
		catch (SingularException) {
			A_tilde_qq = IdentityMatrix(9);
		}
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
	if (Globals::mode == DeformationMode::ClusterBased) {
		// 4. Cluster-Bases deformation
		// calc roofPosition only based on external forces 
		for (u32 index = 0; index < particles.size(); index++) {
			Particle* par = particles[index];
			par->velocityRoof = par->velocity + par->getForces() * timeElapsed;
			par->roofPosition = par->position + par->velocityRoof * timeElapsed;
		}

		// calculate current centre of mass
		core::array<vector3df> currentCentreOfMass, tmpSum;
		core::array<float> tmpMass = 0.0f;
		tmpSum = (0.0f, 0.0f, 0.0f);
		for (int i = 0; i < clusters.size(); i++){
			// this cluster
			vector3df sum(0.0f, 0.0f, 0.0f);
			float mass = 0.0f;
			for (u32 index = 0; index < particles.size(); index++) {
				Particle* particle = particles[index];
				//tmpSum += particle->roofPosition * particle->weight;
				// sigma(m_i * x_i)
				if (which_cluster(particle) == i){
					// if this particle is in this cluster
					sum.X += particle->weight * particle->roofPosition.X;
					sum.Y += particle->weight * particle->roofPosition.Y;
					sum.Z += particle->weight * particle->roofPosition.Z;
					mass += particle->weight;
				}// this particle over
			}// all the particles in this cluster over
			tmpSum.push_back(sum);
			tmpMass.push_back(mass);
			// next cluster
		}// all the clusters over
		// x_cm = sigma(m_i * x_i) / sigma(m_i)
		for (int i = 0; i < clusters.size(); i++){
			vector3df centre(0.0f, 0.0f, 0.0f);
			centre.X += tmpSum[i].X / tmpMass[i];
			centre.Y += tmpSum[i].Y / tmpMass[i];
			centre.Z += tmpSum[i].Z / tmpMass[i];
			currentCentreOfMass.push_back(centre);
		}
	}
	else{
		vector3df tmpSum(0.0f, 0.0f, 0.0f);
		float tmpMass = 0.0f;
		for (u32 index = 0; index < particles.size(); index++) {
			Particle* par = particles[index];
			tmpSum += par->roofPosition * par->weight;
			tmpMass += par->weight;
		}
		currentCentreOfMass = tmpSum / tmpMass;
	}

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
	else if (Globals::mode == DeformationMode::ClusterBased) {
		// 4. Cluster-Bases deformation
		matrix4 A_pq = calculateA_pqMatrix();
		matrix4 R = calculateRotationMatrix(A_pq);

		// calculate goal positions and integrate
		for (u32 i = 0; i < particles.size(); i++) {
			vector3df goalPosition = q[i];
			R.rotateVect(goalPosition);
			goalPosition = goalPosition + currentCentreOfMass;

			if (particles[i]->canMoveTo(goalPosition)) {
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
