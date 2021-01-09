#include <array>
#include <memory>
#include <vector>

#include "ParametrisedTransforms.h"

struct ImMat4 {
	ImVec4 col1, col2, col3, col4;
};
struct ImMat2 {
	ImVec2 col1, col2;
};


typedef int ImGuiControlPointFlags;

enum ImGuiControlPointFlags_
{
	ImGuiControlPointFlags_DrawControlPointMarkers = 1 << 0,
	ImGuiControlPointFlags_DrawParamDerivatives = 1 << 1,
	ImGuiControlPointFlags_DoNotChangeParams = 1 << 2,
	ImGuiControlPointFlags_ApplyParamChangesImmediately = 1 << 3,
	ImGuiControlPointFlags_Circular = 1 << 4,
	ImGuiControlPointFlags_FixedSize = 1 << 5,
	ImGuiControlPointFlags_SizeInPixels = 1 << 6,
};

namespace ImGui {

	// Calls idiomatic to Dear ImGui, introduces a global view class
	// NB this shares transformation state, it is not cleared at the end of the view
	bool BeginControlPointView(const char* name, ImVec2 const& size, ImVec4 const& border_color);
	void EndControlPointView();

	// Typically we want to apply changes at the end of the view, however if we wish to bring some changes forward we can use a
	// deferral slot, changes will be applied when the slot is popped.
	void PushDeferralSlot();
	void PopDeferralSlot();

	ImVec2 GetControlPointViewBoundsMin();
	ImVec2 GetControlPointViewBoundsMax();

	ImVec2 GetLastControlPointPosition();
	float GetLastControlPointRadius();

	bool ControlPoint(ImVec2 const& pos, ImGuiControlPointFlags flags, ImGuiButtonFlags button_flags, float marker_radius, ImVec4 marker_col, float z_order);
	bool ControlPoint(ImVec2 const& pos, float* free_parameter, ImGuiControlPointFlags flags, ImGuiButtonFlags button_flags, float marker_radius, ImVec4 marker_col, float z_order);
	bool ControlPoint(ImVec2 const& pos, float* free_parameter1, float* free_parameter2, ImGuiControlPointFlags flags, ImGuiButtonFlags button_flags, float marker_radius, ImVec4 marker_col, float z_order);
	bool ControlPoint(ImVec2 const& pos, float* free_parameter1, float* free_parameter2, float* free_parameter3, ImGuiControlPointFlags flags, ImGuiButtonFlags button_flags, float marker_radius, ImVec4 marker_col, float z_order);
	
	bool ControlPoint(ImVec4 const& pos, ImGuiControlPointFlags flags, ImGuiButtonFlags button_flags, float marker_radius, ImVec4 marker_col);
	bool ControlPoint(ImVec4 const& pos, float* free_parameter, ImGuiControlPointFlags flags, ImGuiButtonFlags button_flags, float marker_radius, ImVec4 marker_col);
	bool ControlPoint(ImVec4 const& pos, float* free_parameter1, float* free_parameter2, ImGuiControlPointFlags flags, ImGuiButtonFlags button_flags, float marker_radius, ImVec4 marker_col);
	bool ControlPoint(ImVec4 const& pos, float* free_parameter1, float* free_parameter2, float* free_parameter3, ImGuiControlPointFlags flags, ImGuiButtonFlags button_flags, float marker_radius, ImVec4 marker_col);


	void PushMatrix(ImMat4 const& M);  // 3d case
	void PushPerspectiveMatrix(float fov, float aspect_ratio, float z_near, float z_far);
	void PushLookAtMatrix(Vec3<float> const& eye, Vec3<float> const& center, Vec3<float> const& up);
	void PushMatrix(ImMat2 const& M);  // 2d case

	void PushRotationAboutAxis(float* parameter, Vec3<float> const& axis);  // 3d cases
	void PushRotationAboutAxis(float angle, Vec3<float> const& axis);
	void PushRotationAboutX(float* parameter);
	void PushRotationAboutX(float angle);
	void PushRotationAboutY(float* parameter);
	void PushRotationAboutY(float angle);
	void PushRotationAboutZ(float* parameter);
	void PushRotationAboutZ(float angle);
	void PushRotation(float* parameter);  // 2d cases
	void PushRotation(float angle);

	void PushScale(float* parameter);
	void PushScale(float factor);
	void PushScaleX(float* parameter);
	void PushScaleX(float factor);
	void PushScaleY(float* parameter);
	void PushScaleY(float factor);
	void PushScaleZ(float* parameter);
	void PushScaleZ(float factor);

	void PushTranslation(float* parameter, ImVec4 const& v);
	void PushTranslation(float value, ImVec4 const& v);
	void PushTranslation(float* parameter, ImVec2 const& v);
	void PushTranslation(float value, ImVec2 const& v);
	void PushTranslationAlongX(float* parameter);
	void PushTranslationAlongX(float value);
	void PushTranslationAlongY(float* parameter);
	void PushTranslationAlongY(float value);
	void PushTranslationAlongZ(float* parameter);
	void PushTranslationAlongZ(float value);

	void PopTransform();
	void ClearTransforms();


	class DraggableView {
	public:
		bool begin(const char* name, ImVec2 const& size, ImVec4 const& border_color);
		void end();

		void pushDeferralSlot();
		void popDeferralSlot();

		ImVec2 getViewBoundsMin() const;
		ImVec2 getViewBoundsMax() const;
		ImVec2 const& getLastControlPointPosition() const { return m_last_control_point_position; }
		float getLastControlPointRadius() const { return m_last_control_point_radius; }
		float const& getLastParameterStep(size_t ix) const { return m_last_step[ix].change; }

		bool active{ false };

	protected:
		unsigned int m_cp_count{};
		ImVec2 m_initial_cursor_pos{};
		ImVec2 m_final_cursor_pos{};
		ImVec2 m_size{};

		void drawDerivative(ImVec2 const& pos, ImVec2 const& d, ImVec4 const& col, float thickness = 1.0f) const;
		
		ImVec2 getMouseInViewCoords() const;
		ImVec2 viewCoordsToScreenCoords(ImVec2 const& p) const;
		ImVec2 screenCoordsToViewCoords(ImVec2 const& p) const;
		bool createControlPoint(ImVec2 const& marker_pos, ImGuiControlPointFlags const& flags, ImGuiButtonFlags const& button_flags, float marker_width, ImVec4 marker_col, float z_order);

		ImGuiID m_saved_control_point_id{};
		float m_saved_control_point_z_value{};

		bool saveControlPoint(float z_order);

		ImVec2 m_last_control_point_position{};
		float m_last_control_point_radius{};

		struct ParameterChange {
			float* parameter{ nullptr };
			float change{ 0 };
		};
		typedef std::array<ParameterChange, 3> SavedParameterChanges;
		SavedParameterChanges m_last_step{};

		// We use the language of a stack of deferrals, but since there is only one change per frame it
		// is implemented as a single change and two unsigned ints, one for the size of stack and one
		// for the position of the change.
		SavedParameterChanges m_deferred_step{};
		unsigned int m_deferral_stack_size{ 1 };  // Stack is never empty, there is always a slot
		int m_deferred_step_position{ -1 };  // -1 means there is no deferred step in the stack

		void updateParameters(SavedParameterChanges const& changes);
		void updateParameters(SavedParameterChanges const& changes, ImGuiControlPointFlags const& flags);
	};

	class Draggable2DView : public DraggableView {
	public:
		ImVec2 apply(ImVec2 const& pos);
		bool controlPoint(ImVec2 const& pos, 
			ImGuiControlPointFlags flags = 0, ImGuiButtonFlags button_flags = 0, float marker_radius = 0.1f, ImVec4 marker_col = { 1, 1, 1, 1 }, float z_order = 0);
		bool controlPoint(ImVec2 const& pos, float* free_param,
			ImGuiControlPointFlags flags = 0, ImGuiButtonFlags button_flags = 0, float marker_radius = 0.1f, ImVec4 marker_col = { 1, 1, 1, 1 }, float z_order = 0);
		bool controlPoint(ImVec2 const& pos, float* free_param1, float* free_param2,
			ImGuiControlPointFlags flags = 0, ImGuiButtonFlags button_flags = 0, float marker_radius = 0.1f, ImVec4 marker_col = { 1, 1, 1, 1 }, float z_order = 0);
		bool controlPoint(ImVec2 const& pos, float* free_param1, float* free_param2, float* free_param3,
			ImGuiControlPointFlags flags = 0, ImGuiButtonFlags button_flags = 0, float marker_radius = 0.1f, ImVec4 marker_col = { 1, 1, 1, 1 }, float z_order = 0);

		template<unsigned int D>
		bool controlPoint(ImVec2 const& pos, std::array<float*, D> free_parameters, 
			ImGuiControlPointFlags flags, ImGuiButtonFlags button_flags, float marker_radius, ImVec4 marker_col, float z_order);

		void pushMatrix(ImMat2 const& M);
		void pushRotation(float* parameter);
		void pushRotation(float angle);

		void pushScale(float* parameter);
		void pushScale(float factor);
		void pushScaleX(float* parameter);
		void pushScaleX(float factor);
		void pushScaleY(float* parameter);
		void pushScaleY(float factor);

		void pushTranslation(float* parameter, ImVec2 const& v);
		void pushTranslation(float value, ImVec2 const& v);
		void pushTranslationAlongX(float* parameter);
		void pushTranslationAlongX(float value);
		void pushTranslationAlongY(float* parameter);
		void pushTranslationAlongY(float value);

		void popTransform() { m_transforms.simple_transforms.pop_back(); }
		void clearTransforms() { m_transforms.simple_transforms.clear(); }

	private:
		Parametrised2dTransformContext<float, ImVec2, ImMat2> m_transforms{};

		template<unsigned int D>
		void bringTogether(ImVec2 pos, ImVec2 mouse_in_view_coords, std::array<float*, D> const& free_parameters, ImGuiControlPointFlags const& flags);
	};

	class Draggable3DView : public DraggableView {
	public:
		ImVec2 apply(ImVec4 const& pos);
		bool controlPoint(ImVec4 const& pos, 
			ImGuiControlPointFlags flags = 0, ImGuiButtonFlags button_flags = 0, float marker_radius = 0.1f, ImVec4 marker_col = { 1, 1, 1, 1 });
		bool controlPoint(ImVec4 const& pos, float* free_param, 
			ImGuiControlPointFlags flags = 0, ImGuiButtonFlags button_flags = 0, float marker_radius = 0.1f, ImVec4 marker_col = { 1, 1, 1, 1 });
		bool controlPoint(ImVec4 const& pos, float* free_param1, float* free_param2,
			ImGuiControlPointFlags flags = 0, ImGuiButtonFlags button_flags = 0, float marker_radius = 0.1f, ImVec4 marker_col = { 1, 1, 1, 1 });
		bool controlPoint(ImVec4 const& pos, float* free_param1, float* free_param2, float* free_param3,
			ImGuiControlPointFlags flags = 0, ImGuiButtonFlags button_flags = 0, float marker_radius = 0.1f, ImVec4 marker_col = { 1, 1, 1, 1 });

		template<unsigned int D>
		bool controlPoint(ImVec4 const& pos, std::array<float*, D> free_parameters, 
			ImGuiControlPointFlags flags, ImGuiButtonFlags button_flags, float marker_radius, ImVec4 marker_col);

		void pushMatrix(ImMat4 const& M);
		void pushPerspectiveMatrix(float fov, float aspect_ratio, float z_near, float z_far);
		void pushLookAtMatrix(Vec3<float> const& eye, Vec3<float> const& center, Vec3<float> const& up);
		
		void pushRotationAboutAxis(float* parameter, Vec3<float> const& axis);
		void pushRotationAboutAxis(float angle, Vec3<float> const& axis);
		void pushRotationAboutX(float* parameter);
		void pushRotationAboutX(float angle);
		void pushRotationAboutY(float* parameter);
		void pushRotationAboutY(float angle);
		void pushRotationAboutZ(float* parameter);
		void pushRotationAboutZ(float angle);

		void pushScale(float* parameter);
		void pushScale(float factor);
		void pushScaleX(float* parameter);
		void pushScaleX(float factor);
		void pushScaleY(float* parameter);
		void pushScaleY(float factor);
		void pushScaleZ(float* parameter);
		void pushScaleZ(float factor);

		void pushTranslation(float* parameter, ImVec4 const& v);
		void pushTranslation(float value, ImVec4 const& v);
		void pushTranslationAlongX(float* parameter);
		void pushTranslationAlongX(float value);
		void pushTranslationAlongY(float* parameter);
		void pushTranslationAlongY(float value);
		void pushTranslationAlongZ(float* parameter);
		void pushTranslationAlongZ(float value);

		void popTransform() { m_transforms.simple_transforms.pop_back(); }
		void clearTransforms() { m_transforms.simple_transforms.clear(); }


	private:
		ParametrisedProjective3dTransformContext<float, ImVec4, ImMat4> m_transforms{};

		template<unsigned int D>
		void bringTogether(ImVec4 pos, ImVec4 mouse_in_view_coords, std::array<float*, D> const& free_parameters, ImGuiControlPointFlags const& flags);
	};


	void TranslationWidget(ImMat4* M, float axes_size, ImGuiControlPointFlags control_point_flags = 0, bool include_planes = true);
}