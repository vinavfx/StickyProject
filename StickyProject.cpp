/*
 * StickyProject.cpp
 * A Nuke plugin to bake projected UV coordinates onto geometry vertices based on a reference frame
 *
 *  Created on: Jun 22, 2011
 *  Author: Ivan Busquets
 *
 *  Based on the standard UVProject code
 */

#include "DDImage/GeoOp.h"
#include "DDImage/Scene.h"
#include "DDImage/CameraOp.h"
#include "DDImage/Knob.h"
#include "DDImage/Knobs.h"
#include "DDImage/Mesh.h"
#include "DDImage/Channel3D.h"
#include "DDImage/ViewFrustum.h"

#include <assert.h>

using namespace DD::Image;

static const char* const CLASS = "StickyProject";
static const char* const HELP = "Bake a UV projection onto points and vertices for a given reference frame";

// Point/Vertex:
enum { POINT_LIST = 0, VERTEX_LIST };
static const char* const object_types[] = {
		"points", "vertices", 0
};

// Projection type:
enum {
	OFF = 0, PERSPECTIVE, PLANAR, SPHERICAL, CYLINDRICAL
};
static const char* const proj_types[] = {
		"off", "perspective", "planar", "spherical", "cylindrical", 0
};

enum {
	eBoth = 0,
	eFront,
	eBack
};

static const char* const SurfaceNames[] = {
		"both", "front", "back", 0
};

class StickyProject : public GeoOp
{
private:
	int projection;
	double u_scale, v_scale;
	bool u_invert, v_invert;
	int plane;
	Matrix4 xform;
	Matrix4 projectxform;
	Matrix4 _perspxform;

	float inv_u_scale, inv_v_scale;
	int outFrame;

	bool _frustumCulling;   // if true don't affect the uv coordinate outside the frustum
	int _surface;           // surface to project onto
	Vector3 _cam_dir;       // projection vector

protected:
	void _validate(bool for_real)
	{
		// Validate the inputs:
		input0()->validate(for_real);

		// Check if input 1 is connected and get the camera xform from it
		Op* op = Op::input(1,1);
		if (dynamic_cast<CameraOp*>(op)) {
			op->validate(for_real);
			CameraOp* cam = (CameraOp*)op;
			projectxform.translation(0.5f, 0.5f, 0.0f);
			projectxform.scale(0.5f, cam->film_width() / cam->film_height() * 0.5f, 1.0f);
			projectxform *= cam->projection();
			xform = cam->imatrix();
			_cam_dir = cam->matrix().z_axis();
		}
		else if (dynamic_cast<AxisOp*>(op)) {
			op->validate(for_real);
			xform = ((AxisOp*)op)->imatrix();
			projectxform.makeIdentity();
			_cam_dir = -((AxisOp*)op)->matrix().z_axis();
		}
		else {
			xform.makeIdentity();
			projectxform.makeIdentity();
		}

		inv_u_scale = (float)(1.0 / u_scale);
		inv_v_scale = (float)(1.0 / v_scale);

		// Calculate the geometry hashes:
		GeoOp::_validate(for_real);
	}

public:
	static const Description description;
	const char* Class() const
	{
		return CLASS;
	}
	const char* node_help() const
	{
		return HELP;
	}

	StickyProject(Node* node) : GeoOp(node)
	{
		projection = PERSPECTIVE;
		u_scale = v_scale = 1.0;
		u_invert = v_invert = false;
		plane = PLANE_XY;
		_surface = eBoth;
		_frustumCulling = false;
		outFrame = 0;
	}

	int minimum_inputs() const
	{
		return 2;
	}
	int maximum_inputs() const
	{
		return 2;
	}

	const OutputContext& inputContext (int n, int offset, OutputContext & scratch) const {
		scratch = this->outputContext();
		if (offset == 1) {
			scratch.setFrame(outFrame);
		}

		return scratch;

	}

	int split_input(int n) const {
		switch (n) {
		case 0: return 2;
		case 1: return 2;
		}
	}


	bool test_input(int input, Op* op) const
	{
		if (input == 1)
			return dynamic_cast<AxisOp*>(op) != 0;
		return GeoOp::test_input(input, op);
	}

	Op* default_input(int input) const
	{
		if (input == 1)
			return 0;
		return GeoOp::default_input(input);
	}

	const char* input_label(int input, char* buffer) const
	{
		switch (input) {
		case 0:
			return 0;
		case 1:
			return "cam";
		default:
			return 0;
		}
	}

	void append(DD::Image::Hash& hash)
	{
		hash.append(outputContext().frame());
	}


	void knobs(Knob_Callback f)
	{
		GeoOp::knobs(f);
		Enumeration_knob(f, &projection, proj_types, "projection", "projection");
		Obsolete_knob(f, "destination", 0);
		Enumeration_knob(f, &plane, plane_orientation_modes, "plane", "plane");
		Enumeration_knob(f, &_surface, SurfaceNames, "project_on", "project on");
		Bool_knob(f, &_frustumCulling, "frustum_culling", "view frustum culling");
		Tooltip(f, "Project uv onto points and vertices only inside the view frustum.");

		Bool_knob(f, &u_invert, "u_invert", "invert u");
		Bool_knob(f, &v_invert, "v_invert", "invert v");
		Double_knob(f, &u_scale, IRange(0, 2), "u_scale", "u scale");
		Double_knob(f, &v_scale, IRange(0, 2), "v_scale", "v scale");

		Divider(f);
		Int_knob(f, &outFrame, "ref_frame");
		SetFlags(f, Knob::EARLY_STORE);
	}

	int knob_changed(Knob* k)
	{
		knob("plane")->enable(projection > PERSPECTIVE);
		knob("u_scale")->enable(projection > PERSPECTIVE);
		knob("v_scale")->enable(projection > PERSPECTIVE);
		knob("frustum_culling")->enable(projection == PERSPECTIVE);
		return 1;
	}

	/*! Hash up knobs that affect the primitive attributes. */
	void get_geometry_hash()
	{
		// Get all hashes up-to-date
		GeoOp::get_geometry_hash();
		if (projection == OFF)
			return;
		// Hash up knobs that affect the UV attributes
		Hash knob_hash;
		knob_hash.reset();
		// Take transform into account:
		xform.append(knob_hash);
		// Take projection matrix into consideration if perspective mode
		if (projection == PERSPECTIVE)
			projectxform.append(knob_hash);
		// Hash rest of local knobs:
		knob_hash.append(projection);
		knob_hash.append(plane);
		knob_hash.append(u_invert);
		knob_hash.append(v_invert);
		knob_hash.append(u_scale);
		knob_hash.append(v_scale);
		//knob_hash.append(uv_attrib_name);
		knob_hash.append(_surface);
		knob_hash.append(_frustumCulling);
		knob_hash.append(outFrame);
		//knob_hash.append(outputContext().frame());
		// Take point hash into account to force uv upating when positions change
		knob_hash.append(geo_hash[Group_Points]);
		knob_hash.append(geo_hash[Group_Matrix]);

		// Change the point or vertex attributes hash:
		geo_hash[Group_Attributes].append(knob_hash);
		geo_hash[Group_Attributes].append(outputContext().frame());
	}

	inline bool normal_projection_pass_test(Vector3& n) const
	{
		float a = _cam_dir.dot( n );

		if (_surface == eFront){
			if (a < 0) return false;
		}
		else if (_surface == eBack){
			if (a > 0) return false;
		}

		return true;
	}

	// get transformation from geoinfo local to camera
	Matrix4 getLocalToCameraMatrix(const GeoInfo& info)
	{
		Matrix4 m;
		if (info.matrix == Matrix4::identity())
			m = xform;
		else
			m = xform * info.matrix;
		return m;
	}

	void projection_on_polygon(const AttributePtr norm,
			const AttributePtr src,
			Attribute* dst,
			GeoInfo& info,
			DD::Image::GroupType t_group_type,
			DD::Image::GroupType n_group_type )
	{
		const unsigned num_prims = info.primitives();

		// build projection matrix in local camera coordinate
		_perspxform = projectxform * getLocalToCameraMatrix(info);

		// update frustum
		CameraOp* cam = static_cast<CameraOp*>(Op::input(1,1));

		// build camera view frustum in geometry local coordinate
		ViewFrustum frustum;
		if (_frustumCulling){
			frustum.update(*cam, info.matrix.inverse() * cam->matrix());
		}

		// transform all point in world coordinate
		const Vector3* PNTS = info.point_array();

		// vertex index array
		std::vector<unsigned> varray;

		for (unsigned p = 0; p < num_prims; ++p) {

			// get primitive
			const Primitive* prim = info.primitive(p);
			const unsigned num_faces = prim->faces();

			// check all faces in the primitive
			for (unsigned f = 0; f < num_faces; f++){

				// get the number of vertex used by face f
				const unsigned num_verts = prim->face_vertices(f);

				// check polygon against frustum
				if (_frustumCulling){
					// skip polygon
					if (frustum.getVisibility(PNTS, prim, f ) == eNotVisible )
						continue;
				}

				// get all used vertices index
				varray.resize(num_verts);

				if (num_faces == 1){
					for (unsigned i = 0; i < num_verts; i++)
						varray[i] = prim->vertex(i);
				}
				else{
					prim->get_face_vertices(f, &varray[0]);
				}

				// check all vertices in the face
				for (unsigned v = 0; v < num_verts; ++v){
					unsigned vi = varray[v];
					unsigned ti = t_group_type == Group_Points ? vi : prim->vertex_offset() + v;
					unsigned ni = n_group_type == Group_Points ? vi : prim->vertex_offset() + v;

					// avoid to re-test the same vertex again
					if (dst->vector4(ti) != src->vector4(ti))
						continue;

					// get vertex normal
					Vector3 n = norm->normal(ni);

					// rotate normal according world info matrix
					n = info.matrix.vtransform(n);

					// avoid to change the uv coordinate associated with the vertex
					if (normal_projection_pass_test(n))
						project_point(PNTS[vi], dst->vector4(ti));
				}
			}
		}
	}

	void project_point_perspective(int index, GeoInfo& info, GeometryList& out)
	{
		// get normal information
		const AttribContext* N_ref = info.get_attribcontext("N");
		const AttributePtr norm = N_ref ? N_ref->attribute : AttributePtr();

		if(!norm){
			Op::error( "Missing \"N\" channel from geometry");
			return;
		}

		// get the original uv attribute used to restore untouched uv coordinate
		const AttribContext* context = info.get_attribcontext("uv");
		AttributePtr uv_original = context ? context->attribute : AttributePtr();

		if(!uv_original){
			Op::error( "Missing uv channel from geometry" );
			return;
		}

		// we have two possibility:
		// the uv coordinate are stored in Group_Points or in Group_Vertices way
		// the same for normal
		DD::Image::GroupType t_group_type = context->group;  // texture coordinate group type
		DD::Image::GroupType n_group_type = N_ref->group;    // normal group type

		// sanity check
		assert(t_group_type == Group_Points || t_group_type == Group_Vertices);
		assert(n_group_type == Group_Points || n_group_type == Group_Vertices);

		// create a buffer to write on it
		Attribute* uv = out.writable_attribute(index, t_group_type, "uv", VECTOR4_ATTRIB);
		assert(uv);

		// copy all original texture coordinate if available
		if (uv_original){

			// sanity check
			assert(uv->size() == uv_original->size());

			for (unsigned i = 0; i < uv->size(); i++)
				uv->vector4(i) =  uv_original->vector4(i);
		}

		// do the projection on all polygons
		projection_on_polygon( norm, uv_original, uv, info, t_group_type, n_group_type );

	}

	/*! Assign UV attribute to point or vertex attribute list. */
	void geometry_engine(Scene& scene, GeometryList& out)
	{
		GeoOp *staticGeo = static_cast<GeoOp*>(Op::input(0,1));
		GeoOp *movingGeo = static_cast<GeoOp*>(Op::input(0,0));
		movingGeo->get_geometry(scene, out);
		GeometryList staticOut;
		staticGeo->get_geometry(scene, staticOut);
		if (projection == OFF)
			return;

		// the perspective projection need a valid
		if ((projection == PERSPECTIVE) && (Op::input(1)==NULL))
			return;

		// Call the engine on all the caches:
		for (unsigned i = 0; i < out.objects(); i++) {
			GeoInfo& info = staticOut[i];

			if (projection == PERSPECTIVE){
				project_point_perspective( i, info, out );
			}
			else
			{
				// get world transformation
				Matrix4 m = getLocalToCameraMatrix(info);

				// Remove UV vertex attribute, as this takes precedence over a point attribute
				info.delete_group_attribute(Group_Vertices, "uv", VECTOR4_ATTRIB);
				// Create a point attribute
				Attribute* uv = out.writable_attribute(i, Group_Points, "uv", VECTOR4_ATTRIB);
				assert(uv);

				// Project point location and save in UV attribute
				const Vector3* PNTS = info.point_array();
				for (unsigned p = 0; p < info.points(); p++)
					project_point(m.transform(*PNTS++), uv->vector4(p));
			}
		}
	}

	void project_point(const Vector3& in, Vector4& out);
};

#define M_TWOPI M_PI * 2.0
#define DEG2RAD M_PI / 180.0
#define RAD2DEG 180.0 / M_PI

/*! Take the point location and project it back through the camera.
    Where it ends up in the camera aperture is the UV coordinate.
 */
void StickyProject::project_point(const Vector3& in, Vector4& out)
{
	float a, b;
	switch (projection) {
	default:
	case PERSPECTIVE:
		out = _perspxform.transform(in, 1);
		break;
	case PLANAR:
		switch (plane) {
		case PLANE_XY:
			a = in.x;
			b = in.y;
			break;
		case PLANE_YZ:
			a = in.z;
			b = in.y;
			break;
		case PLANE_ZX:
			a = in.x;
			b = in.z;
			break;
		}
		out.set(a * inv_u_scale + 0.5f, b * inv_v_scale + 0.5f, 0, 1);
		break;
		case SPHERICAL: {
			// latitude
			double phi = acos(-in.y);
			// longitude
			double theta = -atan2(-in.x, in.z);
			// Right side
			if (theta <= 0.0)
				theta += M_TWOPI;
			out.set((theta / M_TWOPI) * 0.25f * inv_u_scale, (phi / M_PI - 0.5) * inv_v_scale + 0.5f, 0, 1);
			break;
		}
		case CYLINDRICAL: {
			// longitude
			double theta = -atan2(-in.x, in.z);
			// Right side
			if (theta <= 0.0)
				theta += M_TWOPI;
			out.set((theta / M_TWOPI) * 0.25f * inv_u_scale, in.y * 0.5f * inv_v_scale + 0.5f, 0, 1);
			break;
		}
	}
	if (u_invert)
		out.x = out.w - out.x;
	if (v_invert)
		out.y = out.w - out.y;
}

static Op* build(Node* node)
{
	return new StickyProject(node);
}
const Op::Description StickyProject::description(CLASS, build);

// end of StickyProject.cpp
