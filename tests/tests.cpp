#include "ParametrisedTransforms.h"

#include <glm/glm.hpp>
#include <glm/gtx/rotate_vector.hpp>
#include <glm/gtx/string_cast.hpp>

#include "gtest/gtest.h"

namespace glm {
	std::ostream& operator<<(std::ostream& out, vec2 const& v) { return out << glm::to_string(v); }
	std::ostream& operator<<(std::ostream& out, vec3 const& v) { return out << glm::to_string(v); }
	std::ostream& operator<<(std::ostream& out, vec4 const& v) { return out << glm::to_string(v); }
}

const float PI = 3.14159265359f;

glm::vec3 diff_between_vec4s(glm::vec4 const& a, glm::vec4 const& b)
{
	glm::vec3 a3{ a.x, a.y, a.z };
	a3 /= a.w;
	glm::vec3 b3{ b.x,b.y,b.z };
	b3 /= b.w;
	return a3 - b3;
}

float dist_between_vec4s(glm::vec4 const& a, glm::vec4 const& b)
{
	return glm::length(diff_between_vec4s(a, b));
}



TEST(SimpleTransformTest, Actions) {
	// Some simple tests to ensure that the transformations act on vectors as we expect
	glm::vec3 v{ 1, 2, -1 };
	{  // For translations
		TranslateDirectionTransform<float, glm::vec3> T{ v };
		const float t = 2.0f;
		EXPECT_LE(
			glm::length(T.applyToVector(t, glm::vec3{ 0 }) - t * v),
			1e-15f
		) << "Translation incorrect";
	}

	{  // Rotations
		const float angle = 0.1f;
		RotateXTransform<float, glm::vec3> Rx{};
		RotateYTransform<float, glm::vec3> Ry{};
		RotateZTransform<float, glm::vec3> Rz{};

		EXPECT_LE(
			glm::length(Rx.applyToVector(angle, v) - glm::rotateX(v, angle)),
			1e-10f
		) << "Rotation around x-axis matches glm's implementation";

		EXPECT_LE(
			glm::length(Ry.applyToVector(angle, v) - glm::rotateY(v, angle)),
			1e-10f
		) << "Rotation around y-axis matches glm's implementation";

		EXPECT_LE(
			glm::length(Rz.applyToVector(angle, v) - glm::rotateZ(v, angle)),
			1e-10f
		) << "Rotation around z-axis matches glm's implementation";
	}

	{  // Scalings
		const glm::vec3 v{ 1.2f, 0.2f, -0.3f };
		float const f = 2.0;

		ScaleTransform<float, glm::vec3> S{};
		ScaleAlongAxisTransform<float, glm::vec3, 0> Sx{};
		ScaleAlongAxisTransform<float, glm::vec3, 1> Sy{};
		ScaleAlongAxisTransform<float, glm::vec3, 2> Sz{};

		EXPECT_EQ(
			S.applyToVector(2.0f, v),
			glm::vec3(f * 1.2f, f * 0.2f, f * -0.3f)
		) << "Scale incorrect";

		EXPECT_EQ(
			Sx.applyToVector(2.0f, v), glm::vec3( f * 1.2f, 0.2f, -0.3f )
		) << "Scale in x-axis incorrect";

		EXPECT_EQ(
			Sy.applyToVector(2.0f, v), glm::vec3( 1.2f, f * 0.2f, -0.3f )
		) << "Scale in y-axis incorrect";

		EXPECT_EQ(
			Sz.applyToVector(2.0f, v), glm::vec3( 1.2f, 0.2f, f * -0.3f )
		) << "Scale in z-axis incorrect";
	}
}

float error_in_derivative_approx(BaseSimpleTransform<float, glm::vec3, glm::vec3>* p, float const eps)
{
	const float t = 1.0f;
	const glm::vec3 v{ 0.7f, 0.5f, 0.3f };

	auto residual_error = p->applyToVector(t + eps, v) - p->applyToVector(t, v) - eps * p->dt(t + eps/2, v);

	return glm::length(residual_error);
}

float error_in_second_der_approx(BaseSimpleTransform<float, glm::vec3, glm::vec3>* p, float const eps)
{
	const float t = 1.0f;
	const glm::vec3 v{ 0.7f, 0.5f, 0.3f };

	auto residual_error = p->dt(t + eps, v) - p->dt(t, v) - eps * p->dt2(t + eps/2, v);

	return glm::length(residual_error);
}

TEST(SimpleTransformTest, VectorDerivative) {
	// Performs a test against a numerical approximation

	glm::vec3 v{ 0.1f, 0.2f, -0.3f };

	TranslateDirectionTransform<float, glm::vec3> T{ v };
	RotateAboutAxisTransform<float, glm::vec3> R{ Vec3<float>{v[0], v[1], v[2]} };

	RotateXTransform<float, glm::vec3> Rx{};
	RotateYTransform<float, glm::vec3> Ry{};
	RotateZTransform<float, glm::vec3> Rz{};

	StaticMatrixTransform<float, glm::vec3, glm::vec3, glm::mat3> M{ glm::mat3{2.1f} };
	ScaleTransform<float, glm::vec3> S{};
	ScaleAlongAxisTransform<float, glm::vec3, 0> Sx{};
	ScaleAlongAxisTransform<float, glm::vec3, 1> Sy{};
	ScaleAlongAxisTransform<float, glm::vec3, 2> Sz{};

	float const eps = 1e-4f;
	float const threshold = 20.0f * eps * eps;

	EXPECT_LT(error_in_derivative_approx(&T, eps), threshold) << "Larger error than expected for translation";
	EXPECT_LT(error_in_derivative_approx(&R, eps), threshold) << "Larger error than expected for rotation";
	EXPECT_LT(error_in_derivative_approx(&Rx, eps), threshold) << "Larger error than expected for rotation about x-axis";
	EXPECT_LT(error_in_derivative_approx(&Ry, eps), threshold) << "Larger error than expected for rotation about y-axis";
	EXPECT_LT(error_in_derivative_approx(&Rz, eps), threshold) << "Larger error than expected for rotation about z-axis";
	EXPECT_LT(error_in_derivative_approx(&M, eps), threshold) << "Larger error than expected for static matrix";
	EXPECT_LT(error_in_derivative_approx(&S, eps), threshold) << "Larger error than expected for scale";
	EXPECT_LT(error_in_derivative_approx(&Sx, eps), threshold) << "Larger error than expected for scale along x-axis";
	EXPECT_LT(error_in_derivative_approx(&Sy, eps), threshold) << "Larger error than expected for scale along y-axis";
	EXPECT_LT(error_in_derivative_approx(&Sz, eps), threshold) << "Larger error than expected for scale along z-axis";
}

TEST(SimpleTransformTest, Vector2ndDerivative) {
	// Performs a test against a numerical approximation

	glm::vec3 v{ 0.1f, 0.2f, -0.3f };

	TranslateDirectionTransform<float, glm::vec3> T{ v };
	RotateAboutAxisTransform<float, glm::vec3> R{ Vec3<float>{v[0], v[1], v[2]} };

	RotateXTransform<float, glm::vec3> Rx{};
	RotateYTransform<float, glm::vec3> Ry{};
	RotateZTransform<float, glm::vec3> Rz{};

	StaticAffineTransform<float, glm::vec3, glm::vec3, glm::mat3> M{ glm::mat3{2.1f} };
	ScaleTransform<float, glm::vec3> S{};
	ScaleAlongAxisTransform<float, glm::vec3, 0> Sx{};
	ScaleAlongAxisTransform<float, glm::vec3, 1> Sy{};
	ScaleAlongAxisTransform<float, glm::vec3, 2> Sz{};

	float const eps = 1e-4f;
	float const threshold = 20.0f * eps * eps;

	EXPECT_LT(error_in_second_der_approx(&T, eps), threshold) << "Larger error than expected for translation";
	EXPECT_LT(error_in_second_der_approx(&R, eps), threshold) << "Larger error than expected for rotation";
	EXPECT_LT(error_in_second_der_approx(&Rx, eps), threshold) << "Larger error than expected for rotation about x-axis";
	EXPECT_LT(error_in_second_der_approx(&Ry, eps), threshold) << "Larger error than expected for rotation about y-axis";
	EXPECT_LT(error_in_second_der_approx(&Rz, eps), threshold) << "Larger error than expected for rotation about z-axis";
	EXPECT_LT(error_in_second_der_approx(&M, eps), threshold) << "Larger error than expected for static matrix";
	EXPECT_LT(error_in_second_der_approx(&S, eps), threshold) << "Larger error than expected for scale";
	EXPECT_LT(error_in_second_der_approx(&Sx, eps), threshold) << "Larger error than expected for scale along x-axis";
	EXPECT_LT(error_in_second_der_approx(&Sy, eps), threshold) << "Larger error than expected for scale along y-axis";
	EXPECT_LT(error_in_second_der_approx(&Sz, eps), threshold) << "Larger error than expected for scale along z-axis";
}

TEST(ReparametrisationOfSimpleTransform, TestOfOutput) {

	RotateXTransform<float, glm::vec3> Rx{};
	LinearReparametrisation<RotateXTransform<float, glm::vec3>> RoL{ 2.0f, 1.0f, Rx };

	glm::vec4 v{ 0.1f, 0.2f, -0.3f, 1.2f };
	float t = 1.0f;
	JetDeg2<glm::vec3, float*, 1> J{ {&t, glm::vec3{}} };

	auto RJ = Rx.applyTo2Jet(2 * t + 1, &t, J);
	auto RLJ = Rx.applyTo2Jet(t, &t, J);

	EXPECT_EQ(RJ.position, RLJ.position) << "Position mismatch";
	EXPECT_EQ(RLJ.derivative(0), 2.0f * RJ.derivative(0)) << "Mismatch of derivatives";
	EXPECT_EQ(RLJ.second_derivative(0, 0), 4.0f * RJ.second_derivative(0, 0)) << "Mismatch of second derivatives";

}

TEST(ParametrisedTransform, NumericalTestOfDerivatives) {
	// We numerically test the differentiation carried out by the ParametrisedTransform3dContext class
	// Note that this is not the place to test the individual parameter derivatives, they have their own tests

	float x{ 0.5f }, y{ -0.1f }, z{ 1.1f };

	ParametrisedProjective3dTransformContext<float, glm::vec4, glm::mat4> ctx{};

	ctx.pushRotationAboutX(&x);
	ctx.pushRotationAboutY(&y);
	ctx.pushRotationAboutZ(&z);

	glm::vec4 v{ 0.1f, -2.1f, 1.0f, 1.0f };

	JetDeg1<glm::vec4, float*, 1> initial_1jet{ {&y}, v };
	JetDeg2<glm::vec4, float*, 1> initial_2jet{ initial_1jet };

	// Calculate first and second derivatives
	auto J2 = ctx.applyTo2Jet<1>(initial_2jet);  // Differentiate twice at v with respect to y

	const float eps = 1e-3f;
	const float threshold = 20.0f * eps * eps;
	y += eps;

	auto newJ = ctx.applyTo1Jet<1>(initial_1jet);  // Differentiate at v with respect to y at new value of y

	auto residual = newJ.position - J2.position - eps * J2.derivative(0);
	float error = glm::length(residual);
	EXPECT_LT(error, threshold) << "Large residual error with calculated derivative";

	error = glm::length(newJ.derivative(0) - J2.derivative(0) - eps * J2.second_derivative(0, 0));
	EXPECT_LT(error, threshold) << "Large residual error with calculated second derivative";

}

TEST(ParametrisedTransform, ConsistencyOfDerivatives) {
	// There are multiple ways of computing the same derivatives which should perform the same calculations,
	// make sure they match.

	float x{ 0.5f }, y{ -0.1f }, z{ 1.1f };

	ParametrisedProjective3dTransformContext<float, glm::vec4, glm::mat4> ctx{};

	ctx.pushRotationAboutX(&x);
	ctx.pushRotationAboutY(&y);
	ctx.pushRotationAboutZ(&z);

	glm::vec4 v{ 0.1f, -2.1f, 1.0f, 1.0f };

	JetDeg1<glm::vec4, float*, 1> initial_1jet{ {&y}, v };
	JetDeg2<glm::vec4, float*, 1> initial_2jet{ initial_1jet };

	auto w = ctx.applyToVector(v);
	auto J = ctx.applyTo1Jet<1>(initial_1jet);  // Differentiate at v with respect to y
	auto J2 = ctx.applyTo2Jet<1>(initial_2jet);  // Differentiate twice at v with respect to y

	EXPECT_EQ(w, J.position) << "Base of transformation of 1-jet should be the same as transformation of position.";
	EXPECT_EQ(w, J2.position) << "Base of transformation of 2-jet should be the same as transformation of position.";
	EXPECT_EQ(J.derivative(0), J2.derivative(0)) << "Derivative of transformation of 2-jet should be the same as transformation of derivative.";
}