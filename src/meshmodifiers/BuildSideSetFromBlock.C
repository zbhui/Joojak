
#include "BuildSideSetFromBlock.h"

template<>
InputParameters validParams<BuildSideSetFromBlock>()
{
	InputParameters params = validParams<MeshModifier>();
//	params.addParam<SubdomainName>("from_block", "  ");
	params.addParam<BoundaryName>("new_boundary", " ");
	return params;
}
BuildSideSetFromBlock::BuildSideSetFromBlock(const std::string& name, InputParameters parameters):
    	    MeshModifier(name, parameters)
{
	_boundary_name = "boundary_from_block";
	if(isParamValid("new_boundary"))
		_boundary_name = getParam<BoundaryName>("new_boundary");
}

void BuildSideSetFromBlock::modify()
{
	MeshBase & mesh = _mesh_ptr->getMesh();

	BoundaryInfo & boundary_info = mesh.get_boundary_info();

	MeshBase::const_node_iterator node_it = mesh.active_nodes_begin();
	const MeshBase::const_node_iterator node_end = mesh.active_nodes_end();

	for (; node_it != node_end ; ++node_it)
	{
		const Node* node = *node_it;
		boundary_info.add_node(node, 1000);
	}

	boundary_info.sideset_name(1000) = _boundary_name;
}
