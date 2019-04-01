#include "driver_state.h"
#include <cstring>
using namespace std;
#include <typeinfo>     // for the `std::type_info` type



//
// vertex_shader color
// fragment_shader gouraud
// vertex_data ffffff
// v 0.5 0.25 0.5 1 0 0
// v 0.75 0.25 -0.5 0 1 0
// v 0.75 0.75 -0.5 0 0 1
// render triangle
const data_geometry* ind(driver_state& state, const data_geometry* in[3] );

driver_state::driver_state()
{
}

driver_state::~driver_state()
{
    delete [] image_color;
    delete [] image_depth;
}

// This function should allocate and initialize the arrays that store color and
// depth.  This is not done during the constructor since the width and height
// are not known when this class is constructed.
void initialize_render(driver_state& state, int width, int height)
{
    state.image_width=width;
    state.image_height=height;

    state.image_color= new pixel[width*height];
    state.image_depth= new float[width*height];

    for(int i = 0; i < width*height; i++ ){
        state.image_color[i] = make_pixel(0,0,0);
        state.image_depth[i] = -1;

    }



    std::cout<<"TODO: allocate and initialize state.image_color and state.image_depth."<<std::endl;
}

const data_geometry* Versh(driver_state& state, const data_geometry* in[3] );
void clip_triangle(driver_state& state, const data_geometry out[3],int face);




// This function will be called to render the data that has been stored in this class.
// Valid values of type are:
//   render_type::triangle - Each group of three vertices corresponds to a triangle.
//   render_type::indexed -  Each group of three indices in index_data corresponds
//                           to a triangle.  These numbers are indices into vertex_data.
//   render_type::fan -      The vertices are to be interpreted as a triangle fan.
//   render_type::strip -    The vertices are to be interpreted as a triangle strip.
void render(driver_state& state, render_type type)
{
  int index = 0;
  const data_geometry** tri = new const data_geometry*[3];
  if(type == render_type::triangle ){
    cout << "testing tirangle \n";

  for(int i = 0; i < state.num_vertices; i++){
      data_geometry* creating = new  data_geometry;
        creating->data = new float[state.floats_per_vertex];
      for(int j =0; j < state.floats_per_vertex; j++){ // floats per vertex
          creating->data[j] = state.vertex_data[(i*state.floats_per_vertex)+j];
      }
      tri[index] = creating;
      index++;
      if(index %  3 == 0 ){
        index=0;

        Versh(state,tri);
        tri = new const data_geometry*[3];
      }
  }
}else if(type == render_type::indexed ){
  cout << "testing\n";
  for(int i = 0; i < state.num_triangles*3; i++){
    // cout << "indexing " << state.index_data[i] << endl;

      data_geometry* creating = new  data_geometry;
        creating->data = new float[3];
        // cout << "vertices";
        int wow = 0;
      for(int j =3*state.index_data[i]; j < 3*state.index_data[i] + 3; j++){ // floats per vertex
          creating->data[wow] = state.vertex_data[j];
          wow++;
      }
      tri[index] = creating;
      index++;
      // cout << "the number index " <<  index <<  endl;
      // cout << "the number vertices " << state.num_vertices <<  endl;

      if(index %  3 == 0 ){
        index=0;

        ind(state,tri);
        tri = new const data_geometry*[3];
      }
  }

}else if(type == render_type::strip ){


  for(int i = 0; i < state.num_vertices; i++){
    // cout << "indexing " << state.index_data[i] << endl;

      data_geometry* creating = new  data_geometry;
        creating->data = new float[state.floats_per_vertex];
      for(int j =0; j < state.floats_per_vertex; j++){ // floats per vertex
            creating->data[j] = state.vertex_data[(i*state.floats_per_vertex)+j];
            cout << "the first vertices" <<creating->data[j] <<endl;
      }
      cout << "\n";
      tri[index] = creating;
      index++;
      if(index %  3 == 0 ){
        index=0;
        i = i -2;
        ind(state,tri);
        tri = new const data_geometry*[3];
      }
  }

}else if(type == render_type::fan ){
  cout << "testing\n";
  cout << "num of triangle" << state.num_vertices <<endl;;



  for(int i = 0; i < state.num_vertices; i++){
    // cout << "indexing " << state.index_data[i] << endl;

      data_geometry* creating = new  data_geometry;
        creating->data = new float[state.floats_per_vertex];
      for(int j =0; j < state.floats_per_vertex; j++){ // floats per vertex
          if (index == 0) {
              creating->data[j] = state.vertex_data[j];
          }else{
            creating->data[j] = state.vertex_data[(i*state.floats_per_vertex)+j];

          }
      }
      cout << "\n";
      tri[index] = creating;
      index++;
      if(index %  3 == 0 ){
        index=0;
        i = i -2;
        ind(state,tri);
        tri = new const data_geometry*[3];
      }
  }

}
    std::cout<<"TODO: implement rendering."<<std::endl;
}


const data_geometry* Versh(driver_state& state, const data_geometry* in[3] ){

  data_geometry* out = new  data_geometry[3];

  // initialize data for the output
  out[0].data = new float[state.floats_per_vertex];
  out[1].data = new float[state.floats_per_vertex];
  out[2].data = new float[state.floats_per_vertex];
  // data_geometry* outPorperct = new  data_geometry[3];


  //caling the vertex shader
  state.vertex_shader({in[0]->data},out[0] ,state.uniform_data); // one vertex at the time  ??  or what ?
  state.vertex_shader({in[1]->data},out[1] ,state.uniform_data); // one vertex at the time  ??  or what ?
  state.vertex_shader({in[2]->data},out[2] ,state.uniform_data); // one vertex at the time  ??  or what ?
  // cout << out[0].gl_Position << endl;
  // cout << out[1].gl_Position << endl;
  // cout << out[2].gl_Position << endl;
  const data_geometry** tri = new const data_geometry*[3];
  data_geometry *t1 = new data_geometry;
  data_geometry *t2 = new data_geometry;
  data_geometry *t3 = new data_geometry;
  t1->data = new float[state.floats_per_vertex];
  t2->data = new float[state.floats_per_vertex];
  t3->data = new float[state.floats_per_vertex];

  for(int i=0; i < state.floats_per_vertex; i++){
    t1->data[i] = out[0].data[i];
    t2->data[i] = out[1].data[i];
    t3->data[i] = out[2].data[i];
  }
  //
  t1->gl_Position = out[0].gl_Position;
  t2->gl_Position = out[1].gl_Position;
  t3->gl_Position = out[2].gl_Position;
  tri[0] = t1;
  tri[1] = t2;
  tri[2] = t3;
  clip_triangle(state,tri,0);


}


const data_geometry* ind(driver_state& state, const data_geometry* in[3] ){

  data_geometry* out = new  data_geometry[3];

  // initialize data for the output
  out[0].data = new float[state.floats_per_vertex];
  out[1].data = new float[state.floats_per_vertex];
  out[2].data = new float[state.floats_per_vertex];
  // data_geometry* outPorperct = new  data_geometry[3];


  //caling the vertex shader
  state.vertex_shader({in[0]->data},out[0] ,state.uniform_data); // one vertex at the time  ??  or what ?
  state.vertex_shader({in[1]->data},out[1] ,state.uniform_data); // one vertex at the time  ??  or what ?
  state.vertex_shader({in[2]->data},out[2] ,state.uniform_data); // one vertex at the time  ??  or what ?
  // cout << out[0].gl_Position << endl;
  // cout << out[1].gl_Position << endl;
  // cout << out[2].gl_Position << endl;
  const data_geometry** tri = new const data_geometry*[3];
  data_geometry *t1 = new data_geometry;
  data_geometry *t2 = new data_geometry;
  data_geometry *t3 = new data_geometry;
  t1->data = new float[state.floats_per_vertex];
  t2->data = new float[state.floats_per_vertex];
  t3->data = new float[state.floats_per_vertex];

  for(int i=0; i < state.floats_per_vertex; i++){
    t1->data[i] = out[0].data[i];
    t2->data[i] = out[1].data[i];
    t3->data[i] = out[2].data[i];
  }
  //
  t1->gl_Position = out[0].gl_Position;
  t2->gl_Position = out[1].gl_Position;
  t3->gl_Position = out[2].gl_Position;
  tri[0] = t1;
  tri[1] = t2;
  tri[2] = t3;
  rasterize_triangle(state,tri);


}


// This function clips a triangle (defined by the three vertices in the "in" array).
// It will be called recursively, once for each clipping face (face=0, 1, ..., 5) to
// clip against each of the clipping faces in turn.  When face=6, clip_triangle should
// simply pass the call on to rasterize_triangle.
void clip_triangle(driver_state& state, const data_geometry* out[3],int face)
{
  if(face==6){
    rasterize_triangle(state, out);
    return;
  }




  int  code = 0;
  if(  face == 0){ // checking the cases for clipping
    code |= (out[0]->gl_Position[0] <= abs(out[0]->gl_Position[3]));
    code |= (out[1]->gl_Position[0] <= abs(out[1]->gl_Position[3])) << 1;
    code |= (out[2]->gl_Position[0] <= abs(out[2]->gl_Position[3])) << 2;
  }else if( face  == 1 ){
    code |= (out[0]->gl_Position[0] >= -out[0]->gl_Position[3]);
    code |= (out[1]->gl_Position[0] >= -out[1]->gl_Position[3]) << 1;
    code |= (out[2]->gl_Position[0] >= -out[2]->gl_Position[3]) << 2;
  }else if( face  == 2 ){
    code |= (out[0]->gl_Position[1] <= abs(out[0]->gl_Position[3]));
    code |= (out[1]->gl_Position[1] <= abs(out[1]->gl_Position[3])) << 1;
    code |= (out[2]->gl_Position[1] <= abs(out[2]->gl_Position[3])) << 2;
  }else if( face  == 3 ){
    code |= (out[0]->gl_Position[1] >= -out[0]->gl_Position[3]);
    code |= (out[1]->gl_Position[1] >= -out[1]->gl_Position[3]) << 1;
    code |= (out[2]->gl_Position[1] >= -out[2]->gl_Position[3]) << 2;
  }else if( face  == 4 ){
    code |= (out[0]->gl_Position[2] <= abs(out[0]->gl_Position[3]));
    code |= (out[1]->gl_Position[2] <= abs(out[1]->gl_Position[3])) << 1;
    code |= (out[2]->gl_Position[2] <= abs(out[2]->gl_Position[3])) << 2;
  }else if( face  == 5 ){
    code |= (out[0]->gl_Position[2] >= -out[0]->gl_Position[3]);
    code |= (out[1]->gl_Position[2] >= -out[1]->gl_Position[3]) << 1;
    code |= (out[2]->gl_Position[2] >= -out[2]->gl_Position[3]) << 2;
  }


  float alpha;
  float beta;
  const data_geometry** newTri1 = new  const data_geometry*[3];
  const data_geometry** newTri2 = new const data_geometry*[3];

  float w = 0;
  // //   // initialize data for the output
  data_geometry *t1 = new data_geometry;
  data_geometry *t2 = new data_geometry;
  data_geometry *t3 = new data_geometry;
  t1->data = new float[state.floats_per_vertex];
  t2->data = new float[state.floats_per_vertex];
  t3->data = new float[state.floats_per_vertex];
  data_geometry *t4 = new data_geometry;
  data_geometry *t5 = new data_geometry;
  data_geometry *t6 = new data_geometry;
  t4->data = new float[state.floats_per_vertex];
  t5->data = new float[state.floats_per_vertex];
  t6->data = new float[state.floats_per_vertex];

//  vec4 p;
    switch (code) { //the zero means that  it is outside  1 that it is inside
      case 0b000:
      // cout << "testing " << 0 << endl;;
      return;



      // return;
      case 0b001:
      if (face%2 == 0) {
        w = abs( out[2]->gl_Position[3]);
      }else{
        w = - out[2]->gl_Position[3];
      }
      alpha = (w - out[0]->gl_Position[face/2]) /(out[2]->gl_Position[face/2] - out[0]->gl_Position[face/2] ); // computes it right
      t1->gl_Position  = out[0]->gl_Position;
      t3->gl_Position  = (alpha)*out[2]->gl_Position + (1- alpha)*out[0]->gl_Position;

      // computing the second point and triangle
      if (face%2 == 0) {
        w = abs( out[1]->gl_Position[3]);
      }else{
        w = - out[1]->gl_Position[3];
      }
      beta = (w - out[0]->gl_Position[face/2]) /(out[1]->gl_Position[face/2] - out[0]->gl_Position[face/2] ); // computes it right
      t2->gl_Position  = ((beta))*out[1]->gl_Position +(1- beta)*out[0]->gl_Position;
      for (size_t i = 0; i < state.floats_per_vertex; i++) {
          t1->data[i]= out[0]->data[i];
          t2->data[i]= beta*out[1]->data[i] + (1-beta)*out[0]->data[i];
          t3->data[i]= alpha*out[2]->data[i] + (1-alpha)*out[0]->data[i];
      }

      newTri1[0] = t1;
      newTri1[1] = t2;
      newTri1[2] = t3;

      clip_triangle(state,newTri1,face+1);
      return;





      // cout << "testing " << 1 << endl;

      case 0b010:
      // cout << "testing " << 2 << endl;

      if (face%2 == 0) {
        w = abs( out[0]->gl_Position[3]);
      }else{
        w = - out[0]->gl_Position[3];
      }
      alpha = (w - out[1]->gl_Position[face/2]) /(out[0]->gl_Position[face/2] - out[1]->gl_Position[face/2] ); // computes it right
      t1->gl_Position  = out[1]->gl_Position;
      t2->gl_Position  = (alpha)*out[0]->gl_Position + (1- alpha)*out[1]->gl_Position;
      // computing the  t2->gl_Position  second point and triangle
      if (face%2 == 0) {
        w = abs( out[2]->gl_Position[3]);
      }else{
        w = - out[2]->gl_Position[3];
      }
      beta = (w - out[1]->gl_Position[face/2]) /(out[2]->gl_Position[face/2] - out[1]->gl_Position[face/2] ); // computes it right
      t3->gl_Position  = ((beta))*out[2]->gl_Position +(1- beta)*out[1]->gl_Position;
      for (size_t i = 0; i < state.floats_per_vertex; i++) {
          t1->data[i]= out[1]->data[i];
          t2->data[i]= alpha*out[0]->data[i] + (1-alpha)*out[1]->data[i];
          t3->data[i]= beta*out[2]->data[i] + (1-beta)*out[1]->data[i];
      }
      newTri1[0] = t1;
      newTri1[1] = t2;
      newTri1[2] = t3;

      clip_triangle(state,newTri1,face+1);
      return;





      case 0b011:
      // cout << "testing " << 3 << endl;

      if (face%2 == 0) {
        w = abs( out[2]->gl_Position[3]);
      }else{
        w = - out[2]->gl_Position[3];
      }
      alpha = (w - out[1]->gl_Position[face/2]) /(out[2]->gl_Position[face/2] - out[1]->gl_Position[face/2] ); // computes it right
      t1->gl_Position  = out[0]->gl_Position;
      t2->gl_Position  = out[1]->gl_Position;
      t3->gl_Position  = (alpha)*out[2]->gl_Position + (1- alpha)*out[1]->gl_Position;
      // computing the second point and triangle
      beta = (w - out[0]->gl_Position[face/2]) /(out[2]->gl_Position[face/2] - out[0]->gl_Position[face/2] ); // computes it right
      t4->gl_Position = t3->gl_Position;
      t5->gl_Position  = ((beta))*out[2]->gl_Position +(1- beta)*out[0]->gl_Position;
      t6->gl_Position  = out[0]->gl_Position;

      for (size_t i = 0; i < state.floats_per_vertex; i++) {
          t1->data[i]= out[0]->data[i];
          t2->data[i]= out[1]->data[i];
          t3->data[i]= alpha*out[2]->data[i] + (1-alpha)*out[1]->data[i];

          t4->data[i]= t3->data[i];
          t5->data[i]=  beta*out[2]->data[i] + (1-beta)*out[0]->data[i];
          t6->data[i]= out[0]->data[i];
      }
      newTri1[0] = t1;
      newTri1[1] = t2;
      newTri1[2] = t3;
      newTri2[0] = t4;
      newTri2[1] = t5;
      newTri2[2] = t6;


      clip_triangle(state,newTri1,face+1);
      clip_triangle(state,newTri2,face+1);
      return;



      case 0b100:
      // cout << "testing " << 4 << endl;

      if (face%2 == 0) {
        w = abs( out[0]->gl_Position[3]);
      }else{
        w = - out[0]->gl_Position[3];
      }
      cout << "the w " << w <<endl;
      alpha = (w - out[2]->gl_Position[face/2]) /(out[0]->gl_Position[face/2] - out[2]->gl_Position[face/2] ); // computes it right
      t1->gl_Position  = out[2]->gl_Position;
      t2->gl_Position  = (alpha)*out[0]->gl_Position + (1- alpha)*out[2]->gl_Position;
      // computing the second point and triangle
      // cout << "alfa  " <<   alpha <<endl;
      // cout << "new position 0 - " << out[0]->gl_Position  <<endl;
      // cout << "new position 1 - " << out[1]->gl_Position  <<endl;
      // cout << "new position 2 - " << out[2]->gl_Position  <<endl;

      if (face%2 == 0) {
        w = abs( out[1]->gl_Position[3]);
      }else{
        w = - out[1]->gl_Position[3];
      }
      beta = (w - out[2]->gl_Position[face/2]) /(out[1]->gl_Position[face/2] - out[2]->gl_Position[face/2] ); // computes it right
      t3->gl_Position  = ((beta))*out[1]->gl_Position +(1- beta)*out[2]->gl_Position;
      // cout << "new position  " <<   t2->gl_Position <<endl;
      // cout << "new position  " << t3->gl_Position  <<endl;

      for (size_t i = 0; i < state.floats_per_vertex; i++) {
          t1->data[i]= out[2]->data[i];
          t2->data[i]= alpha*out[0]->data[i] + (1-alpha)*out[2]->data[i];
          t3->data[i]= beta*out[1]->data[i] + (1-beta)*out[2]->data[i];
      }
      newTri1[0] = t1;
      newTri1[1] = t2;
      newTri1[2] = t3;

      rasterize_triangle(state,newTri1);

      // clip_triangle(state,newTri1,face+1);
      return;




      case 0b101:

      // cout << "testing " << 5 << endl;


      if (face%2 == 0) {
        w = abs( out[1]->gl_Position[3]);
      }else{
        w = - out[1]->gl_Position[3];
      }
      alpha = (w - out[0]->gl_Position[face/2]) /(out[1]->gl_Position[face/2] - out[0]->gl_Position[face/2] ); // computes it right
      t1->gl_Position  = out[2]->gl_Position;
      t2->gl_Position  = out[0]->gl_Position;
      t3->gl_Position  = (alpha)*out[1]->gl_Position + (1- alpha)*out[0]->gl_Position;

      // computing the second point and triangle
      beta = (w - out[2]->gl_Position[face/2]) /(out[1]->gl_Position[face/2] - out[2]->gl_Position[face/2] ); // computes it right
      t4->gl_Position = t3->gl_Position;
      t5->gl_Position  = ((beta))*out[1]->gl_Position +(1- beta)*out[2]->gl_Position;

      t6->gl_Position  = out[2]->gl_Position;

      for (size_t i = 0; i < state.floats_per_vertex; i++) {
          t1->data[i]= out[2]->data[i];
          t2->data[i]= out[0]->data[i];
          t3->data[i]= alpha*out[1]->data[i] + (1-alpha)*out[0]->data[i];

          t4->data[i]= t3->data[i];
          t5->data[i]=  beta*out[1]->data[i] + (1-beta)*out[2]->data[i];
          t6->data[i]= out[2]->data[i];
      }
      newTri1[0] = t1;
      newTri1[1] = t2;
      newTri1[2] = t3;
      newTri2[0] = t4;
      newTri2[1] = t5;
      newTri2[2] = t6;

      clip_triangle(state,newTri1,face+1);
      clip_triangle(state,newTri2,face+1);
      return;



      case 0b110:
      // cout << "testing " << 6 << endl;

      if (face%2 == 0) {
        w = abs( out[0]->gl_Position[3]);
      }else{
        w = - out[0]->gl_Position[3];
      }

// (abs(out[0]->gl_Position[3]) - out[2]->gl_Position[face/2]) /(out[0]->gl_Position[face/2] - out[2]->gl_Position[face/2] ); // computes it right
            alpha = (w - out[2]->gl_Position[face/2]) /(out[0]->gl_Position[face/2] - out[2]->gl_Position[face/2] ); // computes it right
            t1->gl_Position  =  out[1]->gl_Position;
            t2->gl_Position  =  out[2]->gl_Position;
            t3->gl_Position  = (alpha)*out[0]->gl_Position + (1- alpha)*out[2]->gl_Position;

// (abs(out[0]->gl_Position[3])*pow(-1,face) - out[1]->gl_Position[face/2]) /(out[0]->gl_Position[face/2] - out[1]->gl_Position[face/2] ); // computes it right
            // computing the second point and triangle
            beta = (w - out[1]->gl_Position[face/2]) /(out[0]->gl_Position[face/2] - out[1]->gl_Position[face/2] ); // computes it right
            t4->gl_Position = t3->gl_Position;
            t5->gl_Position  = ((beta))*out[0]->gl_Position +(1- beta)*out[1]->gl_Position;

            t6->gl_Position  = out[1]->gl_Position;
            for (size_t i = 0; i < state.floats_per_vertex; i++) {
                t1->data[i]= out[1]->data[i];
                t2->data[i]= out[2]->data[i];
                t3->data[i]= alpha*out[0]->data[i] + (1-alpha)*out[2]->data[i];

                t4->data[i]= t3->data[i];
                t5->data[i]=  beta*out[0]->data[i] + (1-beta)*out[1]->data[i];
                t6->data[i]= out[1]->data[i];
            }
            newTri1[0] = t1;
            newTri1[1] = t2;
            newTri1[2] = t3;
            newTri2[0] = t4;
            newTri2[1] = t5;
            newTri2[2] = t6;



            clip_triangle(state,newTri1,face+1);
            clip_triangle(state,newTri2,face+1);
            return;
      case 0b111:
      // cout << "testing " << 0 << endl;

      // cout << "testing " << 7 << endl;
      clip_triangle(state,out, face + 1);
      return;

    }




    std::cout<<"TODO: implement clipping. (The current code passes the triangle through without clipping them.)"<<std::endl;
}
//
//
// // Rasterize the triangle defined by the three vertices in the "in" array.  This
// // function is responsible for rasterization, interpolation of data to
// // fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry *in[3])
{
//   // transform each triangle using the vertex shader
//   // divide by w
//   // transform to pixel coordinates
   data_geometry* out = new  data_geometry[3];

   float wa = in[0]->gl_Position[3];
   float wb = in[1]->gl_Position[3];
   float wc = in[2]->gl_Position[3];
   out[0].gl_Position =  in[0]->gl_Position/ in[0]->gl_Position[3];
   out[1].gl_Position =  in[1]->gl_Position/ in[1]->gl_Position[3];
   out[2].gl_Position=  in[2]->gl_Position/ in[2]->gl_Position[3];


//
    // the transformation
  out[0].gl_Position[0]  = (out[0].gl_Position[0]+1)  * state.image_width/2;
  out[0].gl_Position[1]  = (out[0].gl_Position[1]+1)  * state.image_height/2;
  out[1].gl_Position[0]  = (out[1].gl_Position[0]+1)  * state.image_width/2;
  out[1].gl_Position[1]  = (out[1].gl_Position[1]+1)  * state.image_height/2;
  out[2].gl_Position[0]  = (out[2].gl_Position[0]+1)  * state.image_width/2;
  out[2].gl_Position[1]  = (out[2].gl_Position[1]+1)  * state.image_height/2;

//
   for (size_t y = 0; y <  state.image_height; y++) {

     for(int x = 0; x < state.image_width; x++){
//        /// inter polating checking if it is inside the  triangle.

  // std::cout << "coordinates " <<out[1].gl_Position[0] <<  '\n';

          double abcArea = 0.5 *( (out[1].gl_Position[0]*out[2].gl_Position[1] -  out[2].gl_Position[0]*out[1].gl_Position[1]) -
                   (out[0].gl_Position[0]*out[2].gl_Position[1] -  out[2].gl_Position[0]*out[0].gl_Position[1])+
                   (out[0].gl_Position[0]*out[1].gl_Position[1] -  out[1].gl_Position[0]*out[0].gl_Position[1]) );

          double alpha = 0.5 *( (out[1].gl_Position[0]*out[2].gl_Position[1] -  out[2].gl_Position[0]*out[1].gl_Position[1]) -
                   (x*out[2].gl_Position[1] -  out[2].gl_Position[0]*y)+
                   (x*out[1].gl_Position[1] -  out[1].gl_Position[0]*y) );

          double beta = 0.5 *( (x*out[2].gl_Position[1] -  out[2].gl_Position[0]*y) -
                 (out[0].gl_Position[0]*out[2].gl_Position[1] -  out[2].gl_Position[0]*out[0].gl_Position[1])+
                 (out[0].gl_Position[0]*y -  x*out[0].gl_Position[1]) );

          double delta = 0.5 *( (out[1].gl_Position[0]*y -  x*out[1].gl_Position[1]) -
                 (out[0].gl_Position[0]*y- x*out[0].gl_Position[1]) +
                 (out[0].gl_Position[0]*out[1].gl_Position[1] -  out[1].gl_Position[0]*out[0].gl_Position[1]) );

                 float k  =  ( (alpha/abcArea) /wa  + (beta/abcArea) /wb + (delta/abcArea) /wc );



    float *data_out = new float[state.floats_per_vertex];


        if( (alpha/abcArea > 0) && (delta/abcArea > 0 ) && ( beta/abcArea > 0   ) ){



          data_output* outColor = new  data_output;
          outColor->output_color = vec4(0,0,0,0);
          // different types of interpolations
          for(int i = 0; i  < state.floats_per_vertex; i++){
              if(state.interp_rules[i] == interp_type::flat){ // checking if its a flat
                  data_out[i] = in[0]->data[i];
              } else if(state.interp_rules[i] == interp_type::noperspective){
                  data_out[i] = (alpha/abcArea)*in[0]->data[i] +(beta/abcArea)*in[1]->data[i] + (delta/abcArea)*in[2]->data[i] ;
              }else if(state.interp_rules[i] == interp_type::smooth){
                  data_out[i] = (((alpha/abcArea)/wa)/k)*in[0]->data[i] +(((beta/abcArea)/wb)/k)*in[1]->data[i] + (((delta/abcArea)/wc)/k)*in[2]->data[i] ;
              }
          }

         bool flag =false;
         float  dist =  alpha/abcArea*out[0].gl_Position[2] +beta/abcArea*out[1].gl_Position[2] + delta/abcArea*out[2].gl_Position[2] ;
//          // for z buffering
         if(state.image_depth[ y*state.image_width+x] == -1){
           state.image_depth[ y*state.image_width+x] = dist;
           flag = true;
         }else if( state.image_depth[ y*state.image_width+x] >  dist){
           state.image_depth[ y*state.image_width+x] =dist;
           flag = true;
         }
          state.fragment_shader({data_out},outColor[0] ,state.uniform_data); // one vertex at the time  ??  or what ?

          if(flag){
            // if(beta/abcArea == 0  || alpha/abcArea == 0  || delta/abcArea == 0  ) {
            //   cout << "alpha " <<  alpha/abcArea << endl;
            //   cout << "beta " <<  beta/abcArea << endl;
            //   cout << "delta " <<  delta/abcArea << endl;
            //   cout<< "x y " << x << y << endl;
            //   cout << "colors " <<  outColor[0].output_color[0] << " " << outColor[0].output_color[1]  << " "<< outColor[0].output_color[2] <<endl;
            //   // cout << "multiplaying  " <<  delta/abcArea * 100 << endl;
              state.image_color[ y*state.image_width+x ] = make_pixel( outColor[0].output_color[0]*255,   outColor[0].output_color[1]*255 , outColor[0].output_color[2]*255);


            }
            // state.image_color[ y*state.image_width+x ] = make_pixel( outColor[0].output_color[0]*255,   outColor[0].output_color[1]*255 , outColor[0].output_color[2]*255);
          // }
//
        }

        }
//
   }

    std::cout<<"TODO: implement rasterization"<<std::endl;
}



// this fixed
//
// double abcArea = 0.5 *( (out[1].gl_Position[0]*out[2].gl_Position[1] -  out[2].gl_Position[0]*out[1].gl_Position[1]) -
//          (out[0].gl_Position[0]*out[2].gl_Position[1] -  out[2].gl_Position[0]*out[0].gl_Position[1])+
//          (out[0].gl_Position[0]*out[1].gl_Position[1] -  out[1].gl_Position[0]*out[0].gl_Position[1]) );
//
// double alpha = 0.5 *( (out[1].gl_Position[0]*out[2].gl_Position[1] -  out[2].gl_Position[0]*out[1].gl_Position[1]) -
//          (x*out[2].gl_Position[1] -  out[2].gl_Position[0]*y)+
//          (x*out[1].gl_Position[1] -  out[1].gl_Position[0]*y) );
//
// double beta = 0.5 *( (x*out[2].gl_Position[1] -  out[2].gl_Position[0]*y) -
//        (out[0].gl_Position[0]*out[2].gl_Position[1] -  out[2].gl_Position[0]*out[0].gl_Position[1])+
//        (out[0].gl_Position[0]*y -  x*out[0].gl_Position[1]) );
//
// double delta = 0.5 *( (out[1].gl_Position[0]*y -  x*out[1].gl_Position[1]) -
//        (out[0].gl_Position[0]*y- x*out[0].gl_Position[1]) +
//        (out[0].gl_Position[0]*out[1].gl_Position[1] -  out[1].gl_Position[0]*out[0].gl_Position[1]) );
//

////







///
// double abcArea = 0.5 *( (out[1].gl_Position[0]*out[2].gl_Position[1] -  out[2].gl_Position[0]*out[1].gl_Position[1]) -
//          (out[0].gl_Position[0]*out[2].gl_Position[1] -  out[2].gl_Position[0]*out[0].gl_Position[1])+
//          (out[0].gl_Position[0]*out[1].gl_Position[1] -  out[1].gl_Position[0]*out[0].gl_Position[1]) );
// double alpha = 0.5 *( (out[1].gl_Position[0]*out[2].gl_Position[1] -  out[2].gl_Position[0]*out[1].gl_Position[1]) -
//          (x*out[2].gl_Position[1] -  out[2].gl_Position[0]*y)+
//          (x*out[1].gl_Position[1] -  out[1].gl_Position[0]*y) );
//
// double beta = 0.5 *( (x*out[2].gl_Position[1] -  out[2].gl_Position[0]*y) -
//        (out[0].gl_Position[0]*out[2].gl_Position[1] -  out[2].gl_Position[0]*out[0].gl_Position[1])+
//        (out[0].gl_Position[0]*y -  x*out[0].gl_Position[1]) );
//
// double delta = 0.5 *( (out[1].gl_Position[0]*y -  x*out[1].gl_Position[1]) -
//        (out[0].gl_Position[0]*y- x*out[0].gl_Position[1])+
//        (out[0].gl_Position[0]*out[1].gl_Position[1] -  out[1].gl_Position[0]*out[0].gl_Position[1]) );
//
//        float k  =  ( (alpha/abcArea) /wa  + (beta/abcArea) /wb + (delta/abcArea) /wc );
//
//
