// HyperMesh to Gmsh
// Put it in Altair/12.0/templates/feoutput/
// The components have to be in the same assembly, not in a subassembly, and there should be only one assembly
// The components containing 3D elements, have to contains "fluid" in their name

// Works only with triangle(3 nodes), prism(6 nodes) and tetra(4 nodes)
// To add other types, use 
// http://www.altairhyperworks.com/hwhelp/Altair/hw12.0/help/hwdref/hwdref.aspx?HyperMesh_Reference_Guide.htm
// to find elements(id, ...)
// and http://geuz.org/gmsh/doc/texinfo/gmsh.html#MSH-ASCII-file-format
// to find mshType

// Version of the gmsh file
*text()
    *string("$MeshFormat")
    *end()
    *string("2.2 0 8")
    *end()
    *string("$EndMeshFormat")
    *end()
*output()

*nodes()
  *before()
    *string("$Nodes") *end()
    // Nombre total de noeuds
    *field(integer,[@count(nodes,0,0)],0) *end()
  // id x y z
  *format()
    *field(integer,id,0)
    *string(" ")
    *field(real,x,0)
    *string(" ")
    *field(real,y,0)
    *string(" ")
    *field(real,z,0)
    *end()
  *after()
    *string("$EndNodes")
    *end()
*output()

// tria3
*elements(103,0,"","")
  *before()
    *string("$Elements") *end()
    *field(integer,[@count(elements,0,0)],0) *end()
  // id mshType nbTag idPhEntity idGeoEntity node1-3
  *format()
    *field(integer,id,0) *string(" 2 2 ")
    *field(integer,collector.id,0) *string(" ")
    *field(integer,collector.id,0) *string(" ")
    *field(integer,node1.id,0) *string(" ")
    *field(integer,node2.id,0) *string(" ")
    *field(integer,node3.id,0) *end()
*output()

// penta6
*elements(206,0,"","")
  // id mshType nbTag idPhEntity idGeoEntity node1-6
  *format()
    *field(integer,id,0) *string(" 6 2 ")
    *field(integer,collector.id,0) *string(" ")
    *field(integer,collector.id,0) *string(" ")
    *field(integer,node1.id,0) *string(" ")
    *field(integer,node2.id,0) *string(" ")
    *field(integer,node3.id,0) *string(" ")
    *field(integer,node4.id,0) *string(" ")
    *field(integer,node5.id,0) *string(" ")
    *field(integer,node6.id,0) *end()
*output()

// tetra4
*elements(204,0,"","")
  // id mshType nbTag idPhEntity idGeoEntity node1-4
  *format()
    *field(integer,id,0) *string(" 4 2 ")
    *field(integer,collector.id,0) *string(" ")
    *field(integer,collector.id,0) *string(" ")
    *field(integer,node1.id,0) *string(" ")
    *field(integer,node2.id,0) *string(" ")
    *field(integer,node3.id,0) *string(" ")
    *field(integer,node4.id,0) *end()
  *after()
    *string("$EndElements") *end()
*output()

// Components = Physical Entity
*assemblies()
  *before()
    *string("$PhysicalNames") *end()
  // dim id "name"
  *format()
    *field(integer,numberofcomponents,0) *end()
    *counterset(counter1,0)
    *loopif([counter1 != numberofcomponents])
      *pointerset(pointer1,components,counter1)
      // If the component name contains fluid, it is 3D
      *if([@stringinstring(pointer1.component.name,"fluid",1) == 1])
        *string("3 ")
      // else it is 2D
      *else()
        *string("2 ")
      *endif()
      *field(integer,pointer1.pointervalue,0)
      *string(" ")
      *field(quotedstring,pointer1.component.name,0)
      *end()
      *counterinc(counter1)
    *endloop()
    *after()
      *string("$EndPhysicalNames") *end()
*output()

