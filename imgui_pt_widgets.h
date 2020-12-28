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

namespace ImGui {

	typedef int ControlPointFlags;

	enum ControlPointFlags_
	{
		ControlPointFlags_DrawControlPointMarkers = 1 << 0,
		ControlPointFlags_DrawParameterDerivatives = 1 << 1,
		ControlPointFlags_DoNotChangeParameters = 1 << 2,
		ControlPointFlags_DeferParameterChangeToEnd = 1 << 3,
	};

	class DraggableView {
	public:
		bool begin(const char* name, ImVec2 const& size, ImVec4 const& border_color);
		void end();

		ImVec2 const& getLastControlPointPosition() const { return m_last_control_point_position; }
		float const& getLastParameterStep(size_t ix) const { return m_latest_step[ix].change; }

	protected:
		unsigned int m_cp_count{};
		ImVec2 m_initial_cursor_pos{};
		ImVec2 m_final_cursor_pos{};
		ImVec2 m_size{};

		void drawDerivative(ImVec2 const& pos, ImVec2 const& d, ImVec4 const& col, float thickness = 1.0f) const;
		
		ImVec2 getMouseInViewCoords() const;
		ImVec2 viewCoordsToScreenCoords(ImVec2 const& p) const;
		ImVec2 screenCoordsToViewCoords(ImVec2 const& p) const;
		bool createControlPoint(ImVec2 const& marker_pos, ControlPointFlags const& flags, ImGuiButtonFlags const& button_flags, float marker_width, ImVec4 marker_col, float z_order);

		ImGuiID m_saved_control_point_id{};
		float m_saved_control_point_z_value{};

		bool saveControlPoint(float z_order);

		ImVec2 m_last_control_point_position{};

		struct ParameterChange {
			float* parameter{ nullptr };
			float change{ 0 };
		};
		typedef std::array<ParameterChange, 2> SavedParameterChanges;
		SavedParameterChanges m_latest_step{};
		SavedParameterChanges m_deferred_step{};

		void updateParameters(SavedParameterChanges const& changes, ControlPointFlags const& flags);
	};

	class Draggable2DView : public DraggableView {
	public:
		ImVec2 apply(ImVec2 const& pos);
		bool controlPoint(ImVec2 const& pos, 
			ControlPointFlags flags = 0, ImGuiButtonFlags button_flags = 0, float marker_radius = 0.1f, ImVec4 marker_col = { 1, 1, 1, 1 }, float z_order = 0);
		bool controlPoint(ImVec2 const& pos, float* free_param,
			ControlPointFlags flags = 0, ImGuiButtonFlags button_flags = 0, float marker_radius = 0.1f, ImVec4 marker_col = { 1, 1, 1, 1 }, float z_order = 0);
		bool controlPoint(ImVec2 const& pos, float* free_param1, float* free_param2,
			ControlPointFlags flags = 0, ImGuiButtonFlags button_flags = 0, float marker_radius = 0.1f, ImVec4 marker_col = { 1, 1, 1, 1 }, float z_order = 0);

		template<unsigned int D>
		bool controlPoint(ImVec2 const& pos, std::array<float*, D> free_parameters, 
			ControlPointFlags flags, ImGuiButtonFlags button_flags, float marker_radius, ImVec4 marker_col, float z_order);

		// These are declared here, then defined elsewhere so that the vector functions can be contained in their
		// own translation unit.
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
		void bringTogether(ImVec2 pos, ImVec2 mouse_in_view_coords, std::array<float*, D> const& free_parameters, ControlPointFlags const& flags);
	};

	class Draggable3DView : public DraggableView {
	public:
		ImVec2 apply(ImVec4 const& pos);
		bool controlPoint(ImVec4 const& pos, 
			ControlPointFlags flags = 0, ImGuiButtonFlags button_flags = 0, float marker_radius = 0.1f, ImVec4 marker_col = { 1, 1, 1, 1 });
		bool controlPoint(ImVec4 const& pos, float* free_param, 
			ControlPointFlags flags = 0, ImGuiButtonFlags button_flags = 0, float marker_radius = 0.1f, ImVec4 marker_col = { 1, 1, 1, 1 });
		bool controlPoint(ImVec4 const& pos, float* free_param1, float* free_param2, 
			ControlPointFlags flags = 0, ImGuiButtonFlags button_flags = 0, float marker_radius = 0.1f, ImVec4 marker_col = { 1, 1, 1, 1 });

		template<unsigned int D>
		bool controlPoint(ImVec4 const& pos, std::array<float*, D> free_parameters, 
			ControlPointFlags flags, ImGuiButtonFlags button_flags, float marker_radius, ImVec4 marker_col);

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
		void bringTogether(ImVec4 pos, ImVec4 mouse_in_view_coords, std::array<float*, D> const& free_parameters, ControlPointFlags const& flags);
	};
}