#include "ModelObject.h"

using namespace Model;



//CONSTRUCTORS
ModelObject::ModelObject(ObjectType::ObjectType type) : scene::ISceneNode(Globals::sceneManager->getRootSceneNode(), Globals::sceneManager, (type | ObjectType::ModelObjectType))
{
	// getID() 
	// Get the id of the scene node. This id can be used to identify the node.
	this->type = (ObjectType::ObjectType) getID();
	this->isSelected = false;
	
}

bool ModelObject::is(ObjectType::ObjectType type) {
	return (type & this->type) > 0;
}


// Determine whether we should draw this node, and the children nodes
void ModelObject::OnPreRender()
{
	if (IsVisible) {
		// Draw this node
		// Registers a node for rendering it at a specific time.
		SceneManager->registerNodeForRendering(this);


		// Draw the children nodes
		// Called after the scene's light list has been built, but before rendering has begun.
		ISceneNode::OnPreRender();

	}

}

// Draws this node to the screen
void ModelObject::render()
{


}

// aabbox3d
// Axis aligned bounding box in 3d dimensional space.
// getBoundingBox()
// Get the axis aligned, not transformed bounding box of this node.
const core::aabbox3d<f32>& ModelObject::getBoundingBox() const
{
	// The bounding box of this mesh.
	return Box;
}


