#include "pt_demos.h"

#include "imgui.h"
#include "imgui_pt_widgets.h"

inline float min(float a, float b) { return a < b ? a : b; }
inline float max(float a, float b) { return a > b ? a : b; }
inline float clamp(float x, float a, float b) { return min(max(a, x), b); }

struct TreeParameters {
    float main_length;
    float left_angle;
    float left_scale_factor;
    float right_angle;
    float right_scale_factor;
};

void create_tree(ImGui::Draggable2DView& view, TreeParameters* params, int depth, int LR, ImGuiControlPointFlags flags)
{
    if (depth <= 0)
        return;

    // Branch
    ImVec2 start_point = view.apply({ 0, 0 });
    view.pushTranslationAlongY(&(params->main_length));
    if (LR == 0)
        view.controlPoint({ 0.0f, 0.0f }, &(params->main_length), flags, 0, 0.1f, { 0, 0, 1, 1 }, (float)depth);
    else if (LR == 1) // Right hand branch
        view.controlPoint({ 0.0f, 0.0f }, &(params->right_angle), &(params->right_scale_factor), flags, 0, 0.1f, { 0, 1, 0, 1 }, (float)depth);
    else if (LR == -1) // Left hand branch
        view.controlPoint({ 0.0f, 0.0f }, &(params->left_angle), &(params->left_scale_factor), flags, 0, 0.1f, { 1, 0, 0, 1 }, (float)depth);
    ImVec2 end_point = view.getLastControlPointPosition();

    // Draw the branch
    ImGui::GetWindowDrawList()->AddLine(start_point, end_point, ImGui::GetColorU32({ 1.0, 1.0, 0.0, 1.0 }), 1.0f);

    // Left branch
    view.pushRotation(&(params->left_angle));
    view.pushScale(&(params->left_scale_factor));
    create_tree(view, params, depth - 1, -1, flags);
    view.popTransform();
    view.popTransform();

    // Right branch
    view.pushRotation(&(params->right_angle));
    view.pushScale(&(params->right_scale_factor));
    ImVec2 right_branch = view.getLastControlPointPosition();
    create_tree(view, params, depth - 1, 1, flags);
    view.popTransform();
    view.popTransform();
    view.popTransform();
}

void ShowTreeDemo(bool* p_open)
{
    ImGui::Begin("Fractal Tree - using 2D transforms", p_open);

    static int depth = 5;
    static TreeParameters parameters{ 0.5f, -0.1f, 0.75f, 0.6f, 0.4f };
    static bool draw_markers = false;

    constexpr float factor_lower_bound = 0.2f;
    constexpr float factor_upper_bound = 1.0f;
    constexpr float angle_lower_bound = -3.1416f / 2;
    constexpr float angle_upper_bound = 3.1416f / 2;

    ImGui::SliderFloat("left branch factor", &parameters.left_scale_factor, factor_lower_bound, factor_upper_bound);
    ImGui::SliderFloat("left branch angle", &parameters.left_angle, angle_lower_bound, angle_upper_bound);
    ImGui::SliderFloat("right branch factor", &parameters.right_scale_factor, factor_lower_bound, factor_upper_bound);
    ImGui::SliderFloat("right branch angle", &parameters.right_angle, angle_lower_bound, angle_upper_bound);
    ImGui::SliderFloat("stem length", &parameters.main_length, 0.1f, 1.5f);
    ImGui::SliderInt("tree depth", &depth, 1, 9);
    ImGui::Checkbox("draw control points", &draw_markers);

    static ImGui::Draggable2DView view{};
    if (view.begin("tree_view", ImVec2{ 300, 400 }, ImVec4(0, 0, 0, 1)))
    {
        ImGuiControlPointFlags flags = 0;
        if (draw_markers)
            flags |= ImGuiControlPointFlags_DrawControlPointMarkers;

        view.pushScale(-1.0f);
        view.pushTranslationAlongY(-0.8f);
        create_tree(view, &parameters, depth, 0, flags);
        view.popTransform();
        view.popTransform();

        view.end();  // Any changes to parameters happen in this call, because we deferred the changes

        // Enforce bounds on parameters
        parameters.left_scale_factor = clamp(parameters.left_scale_factor, factor_lower_bound, factor_upper_bound);
        parameters.right_scale_factor = clamp(parameters.right_scale_factor, factor_lower_bound, factor_upper_bound);
        parameters.left_angle = clamp(parameters.left_angle, angle_lower_bound, angle_upper_bound);
        parameters.right_angle = clamp(parameters.right_angle, angle_lower_bound, angle_upper_bound);
    }
    ImGui::End();
}

void ShowSpiralDemo(bool* p_open)
{
    ImGui::Begin("Spiral Demo", p_open);

    static float y_angle = 0.3f;
    static float spiral_rotation = 0.2f;
    static float spiral_height = 0.01f;
    static float spiral_scale = 0.9f;
    static int spiral_length = 20;

    ImGui::SliderFloat("spiral rotation", &spiral_rotation, -3.1416f, 3.1416f);
    ImGui::SliderFloat("spiral height", &spiral_height, -0.1f, 0.1f);
    ImGui::SliderFloat("spiral scale", &spiral_scale, 0.7f, 1.1f);
    ImGui::SliderInt("number of points", &spiral_length, 2, 100);

    ImGui::Draggable3DView view{};

    if (view.begin("some_label_for_id", ImVec2(300, 300), ImVec4(0, 1, 1, 1)))
    {
        view.pushPerspectiveMatrix(45.0f / 180.f * 3.14159f, 1.0f, 0.01f, 10.0f);
        view.pushLookAtMatrix(
            { 2.5f, 0.0f, 1.0f }, { 0.0f, 0.0f, 0.0f }, { 0.0f, 0.0f, -1.0f }
        );
        view.pushRotationAboutY(&y_angle);

        auto draw_list = ImGui::GetWindowDrawList();
        view.controlPoint({ 1.0, 0.0, 0.0, 1.0 }, 0, 0, 0.2f, ImVec4{ 0, 1, 1, 1 });
        ImVec2 prev_point = view.getLastControlPointPosition();
        for (int i = 0; i < spiral_length; ++i)
        {
            view.pushRotationAboutZ(&spiral_rotation);
            view.pushTranslationAlongZ(&spiral_height);
            view.pushScaleX(&spiral_scale);
            view.pushScaleY(&spiral_scale);
            view.controlPoint({ 1.0, 0.0, 0.0, 1.0 }, &spiral_rotation, &spiral_scale, 0, 0, 0.2f, ImVec4(1, 0, 0, 1));
            ImVec2 current_point = view.getLastControlPointPosition();
            draw_list->PathLineTo(prev_point);
            draw_list->PathLineTo(current_point);
            draw_list->PathStroke(ImGui::GetColorU32({ 1.0, 0.0, 0.0, 1.0 }), false, 2.0f);
            prev_point = current_point;
        }
        draw_list->PathStroke(ImGui::GetColorU32({ 1.0, 0.0, 0.0, 1.0 }), false);
        view.pushTranslationAlongZ(0.1f);
        view.pushScaleX(0.0f);
        view.pushScaleY(0.0f);

        ImGuiControlPointFlags flags = ImGuiControlPointFlags_DrawControlPointMarkers;
        view.controlPoint({ 1.0, 0.0, 0.0, 1.0 }, &spiral_height, flags, 0, 0.2f, ImVec4(0, 1, 0, 1));

        view.clearTransforms();
        view.end();
    }

    bool selected = false;
    ImGui::Selectable("Some text to make sure we haven't broken anything", selected);

    ImGui::End();
}

void ShowSpiralDemo2(bool* p_open)
{
    ImGui::Begin("Spiral Demo", p_open);

    static float y_angle = 0.3f;
    static float spiral_rotation = 0.2f;
    static float spiral_height = 0.01f;
    static float spiral_scale = 0.9f;
    static int spiral_length = 20;

    ImGui::SliderFloat("spiral rotation", &spiral_rotation, -3.1416f, 3.1416f);
    ImGui::SliderFloat("spiral height", &spiral_height, -0.1f, 0.1f);
    ImGui::SliderFloat("spiral scale", &spiral_scale, 0.7f, 1.1f);
    ImGui::SliderInt("number of points", &spiral_length, 2, 100);

    if (ImGui::BeginControlPointView("some_label_for_id", ImVec2(300, 300), ImVec4(0, 1, 1, 1)))
    {
        ImGui::PushPerspectiveMatrix(45.0f / 180.f * 3.14159f, 1.0f, 0.01f, 10.0f);
        ImGui::PushLookAtMatrix(
            { 2.5f, 0.0f, 1.0f }, { 0.0f, 0.0f, 0.0f }, { 0.0f, 0.0f, -1.0f }
        );
        ImGui::PushRotationAboutY(&y_angle);

        auto draw_list = ImGui::GetWindowDrawList();
        ImGui::ControlPoint({ 1.0, 0.0, 0.0, 1.0 }, 0, 0, 0.05f, ImVec4{ 0, 1, 1, 1 });
        ImVec2 prev_point = ImGui::GetLastControlPointPosition();
        for (int i = 0; i < spiral_length; ++i)
        {
            ImGui::PushRotationAboutZ(&spiral_rotation);
            ImGui::PushTranslationAlongZ(&spiral_height);
            ImGui::PushScaleX(&spiral_scale);
            ImGui::PushScaleY(&spiral_scale);
            ImGui::ControlPoint({ 1.0, 0.0, 0.0, 1.0 }, &spiral_rotation, &spiral_scale, 0, 0, 0.05f, ImVec4(1, 0, 0, 1));
            ImVec2 current_point = ImGui::GetLastControlPointPosition();
            draw_list->PathLineTo(prev_point);
            draw_list->PathLineTo(current_point);
            draw_list->PathStroke(ImGui::GetColorU32({ 1.0, 0.0, 0.0, 1.0 }), false, 2.0f);
            prev_point = current_point;
        }
        draw_list->PathStroke(ImGui::GetColorU32({ 1.0, 0.0, 0.0, 1.0 }), false);
        ImGui::PushTranslationAlongZ(0.1f);
        ImGui::PushScaleX(0.0f);
        ImGui::PushScaleY(0.0f);

        ImGuiControlPointFlags flags = ImGuiControlPointFlags_DrawControlPointMarkers;
        ImGui::ControlPoint({ 1.0, 0.0, 0.0, 1.0 }, &spiral_height, flags, 0, 0.05f, ImVec4(0, 1, 0, 1));

        ImGui::ClearTransforms();
        ImGui::EndControlPointView();
    }

    bool selected = false;
    ImGui::Selectable("Some text to make sure we haven't broken anything", selected);

    ImGui::End();
}

void ShowArmDemo(bool* p_open)
{
    ImGui::Begin("Arm Demo", p_open);

    static float base_rotation = 0.0f;
    static float base_angle = 0.1f;
    static float elbow_angle = 0.2f;

    ImGuiControlPointFlags flags = ImGuiControlPointFlags_DrawControlPointMarkers;
    flags |= ImGuiControlPointFlags_DrawParamDerivatives;

    ImGui::BeginControlPointView("arm", { 400, 400 }, { 1, 1, 1, 1 });
    ImGui::PushPerspectiveMatrix(45.0f / 180.0f * 3.14159f, 1.0f, 0.01f, 10.0f);
    ImGui::PushLookAtMatrix({ 0.0f, 0.5f, -2.0f }, {}, { 0, -1, 0 });
    ImGui::PushTranslationAlongY(-0.5f);
    ImGui::PushRotationAboutY(&base_rotation);
    ImGui::ControlPoint({ 0.2f, 0, 0, 1.0f }, &base_rotation, flags, 0, 0.025f, { 1, 0, 1, 1 });
    ImGui::ControlPoint({ 0.0f, 0, 0, 1.0f }, flags, 0, 0.025f, { 1, 1, 1, 1 });
    ImGui::PushRotationAboutZ(&base_angle);
    ImGui::PushTranslationAlongY(0.5f);
    ImGui::ControlPoint({ 0.0f, -0.1f, 0, 1.0f }, &base_angle, flags, 0, 0.025f, { 1, 1, 1, 1 });
    ImGui::ControlPoint({ 0.0f, -0.2f, 0, 1.0f }, &base_angle, flags, 0, 0.025f, { 1, 1, 1, 1 });
    ImGui::ControlPoint({ 0.0f, -0.3f, 0, 1.0f }, &base_angle, flags, 0, 0.025f, { 1, 1, 1, 1 });
    ImGui::ControlPoint({ 0.0f, -0.4f, 0, 1.0f }, &base_angle, flags, 0, 0.025f, { 1, 1, 1, 1 });
    ImGui::ControlPoint({ 0.0f, 0, 0, 1.0f }, &base_angle, flags, 0, 0.1f, { 1, 1, 1, 1 });
    ImGui::PushRotationAboutZ(&elbow_angle);
    ImGui::PushTranslationAlongY(0.4f);
    ImGui::ControlPoint({ 0.0f, -0.1f, 0, 1.0f }, &elbow_angle, flags, 0, 0.025f, { 1, 1, 0, 1 });
    ImGui::ControlPoint({ 0.0f, -0.2f, 0, 1.0f }, &elbow_angle, flags, 0, 0.025f, { 1, 1, 0, 1 });
    ImGui::ControlPoint({ 0.0f, -0.3f, 0, 1.0f }, &elbow_angle, flags, 0, 0.025f, { 1, 1, 0, 1 });
    ImGui::ControlPoint({ 0.0f, 0, 0, 1.0f }, &elbow_angle, &base_angle, flags, 0, 0.0375f, { 1, 0, 0, 1 });

    ImGui::ClearTransforms();
    ImGui::EndControlPointView();

    elbow_angle = clamp(elbow_angle, 0.3f, 2.8f);
    base_angle = clamp(base_angle, 0.1f, 1.4f);

    ImGui::End();
}

void ShowTranslationWidgetDemo(bool* p_open)
{
    static ImMat4 model_matrix{
        {1, 0, 0, 0},
        {0, 1, 0, 0},
        {0, 0, 1, 0},
        {0, 0, 0, 1}
    };
    
    ImGui::Begin("Translation Widget", p_open);

    if (ImGui::BeginControlPointView("translation", { 300, 300 }, { 1, 1, 1, 1 }))
    {
        ImGui::PushPerspectiveMatrix(45.0f / 180.0f * 3.14159f, 1.0f, 0.01f, 10.0f);
        ImGui::PushLookAtMatrix({ 0.0f, 0.25f, -1.0f }, {}, { 0, -1, 0 });

        // Draw our control points
        ImGuiControlPointFlags flags = ImGuiControlPointFlags_DrawControlPointMarkers;
        ImGui::TranslationWidget(&model_matrix, 0.1f, flags, true);

        ImGui::PopTransform();
        ImGui::PopTransform();
        ImGui::EndControlPointView();
    }
    if (ImGui::Button("Reset Position")) {
        model_matrix.col4 = { 0, 0, 0, 1 };
    }
    ImGui::End();

}

void ShowCameraDemo(bool* p_open)
{
    ImGui::Begin("Camera Demo", p_open);

    static float cube_rotation = 0.0f;
    static float camera_vertical_angle = 0.0f;
    static float camera_distance = -2.5f;
    static bool show_parameter_derivatives = true;

    ImGui::SliderFloat("cube rotation", &cube_rotation, -10.0f, 10.0f);
    ImGui::SliderFloat("vertical angle", &camera_vertical_angle, -3.1416f / 4, 3.1416f / 4);
    ImGui::SliderFloat("camera displacement", &camera_distance, -5.0f, -1.0f);
    ImGui::Checkbox("show parameter derivatives", &show_parameter_derivatives);

    ImGui::Draggable3DView view{};

    if (view.begin("camera", ImVec2(400, 300), ImVec4(0, 1, 1, 1)))
    {
        ImGuiControlPointFlags flags = show_parameter_derivatives ? ImGuiControlPointFlags_DrawParamDerivatives : 0;
        flags |= ImGuiControlPointFlags_DrawControlPointMarkers;

        view.pushPerspectiveMatrix(45.0f / 180.f * 3.14159f, 4.0f / 3, 0.01f, 10.0f);  // projection matrix

        view.pushTranslationAlongZ(&camera_distance);  // move camera back
        view.pushRotationAboutX(&camera_vertical_angle);
        view.pushRotationAboutY(&cube_rotation);  // Rotation of model along long axis (y is up)
        view.pushLookAtMatrix(  // Model matrix, centering cube at origin and standing it on its diagonal axis
            { 0.5f, 0.5f, 0.5f }, { 1.5f, 0.0f, 0.0f }, { -1.0f, -1.0f, -1.0f }
        );

        // Control points at vertices of cube
        view.controlPoint({ 1.0f, 1.0f, 0.0f, 1.0f }, &cube_rotation, flags, 0, 0.2f, { 1, 0, 0, 1 });
        view.controlPoint({ 1.0f, 0.0f, 1.0f, 1.0f }, &cube_rotation, flags, 0, 0.2f, { 1, 0, 0, 1 });
        view.controlPoint({ 0.0f, 1.0f, 1.0f, 1.0f }, &cube_rotation, flags, 0, 0.2f, { 1, 0, 0, 1 });
        view.controlPoint({ 1.0f, 0.0f, 0.0f, 1.0f }, &camera_vertical_angle, flags, 0, 0.2f, { 1, 0, 0, 1 });
        view.controlPoint({ 0.0f, 1.0f, 0.0f, 1.0f }, &camera_vertical_angle, flags, 0, 0.2f, { 1, 0, 0, 1 });
        view.controlPoint({ 0.0f, 0.0f, 1.0f, 1.0f }, &camera_vertical_angle, flags, 0, 0.2f, { 1, 0, 0, 1 });
        view.controlPoint({ 1.0f, 1.0f, 1.0f, 1.0f }, &camera_distance, flags, 0, 0.2f, { 1, 1, 1, 1 });
        view.controlPoint({ 0.0f, 0.0f, 0.0f, 1.0f }, flags, 0, 0.2f, { 1, 1, 0, 1 });

        view.end();
    }
    ImGui::End();
}