#include <raylib.h>
#include <math.h>
#include <stdio.h>

typedef enum { RK1, RK2, RK3, RK4 } Method;
static Method currentMethod = RK4;

static const float G = 1.0f;
static const float TIME_STEP = 5.0f;

typedef struct {
    Vector2 position;
    Vector2 velocity;
    float mass;
} PointMass;

Vector2 ComputeGravitationalForce(const PointMass* a, const PointMass* b)
{
    const Vector2 direction = (Vector2){ b->position.x - a->position.x, b->position.y - a->position.y };
    float distance = sqrtf(direction.x * direction.x + direction.y * direction.y);
    if (distance < 1.0f) distance = 1.0f; // Avoid division by zero
    const float forceMagnitude = (G * a->mass * b->mass) / (distance * distance);
    return (Vector2){ direction.x / distance * forceMagnitude, direction.y / distance * forceMagnitude };
}

Vector2 ComputeAcceleration(const PointMass* body1, const PointMass* body2)
{
    const Vector2 force = ComputeGravitationalForce(body1, body2);
    return (Vector2){ force.x / body1->mass, force.y / body1->mass };
}

void UpdateRK1(PointMass* body1, PointMass* body2, const float dt)
{
    const Vector2 forceOn1 = ComputeGravitationalForce(body1, body2);
    const Vector2 forceOn2 = (Vector2){ -forceOn1.x, -forceOn1.y };

    const Vector2 acceleration1 = (Vector2){ forceOn1.x / body1->mass, forceOn1.y / body1->mass };
    const Vector2 acceleration2 = (Vector2){ forceOn2.x / body2->mass, forceOn2.y / body2->mass };

    body1->velocity.x += acceleration1.x * dt;
    body1->velocity.y += acceleration1.y * dt;
    body2->velocity.x += acceleration2.x * dt;
    body2->velocity.y += acceleration2.y * dt;

    body1->position.x += body1->velocity.x * dt;
    body1->position.y += body1->velocity.y * dt;
    body2->position.x += body2->velocity.x * dt;
    body2->position.y += body2->velocity.y * dt;
}

void UpdateRK2_Midpoint(PointMass* body1, PointMass* body2, const float dt)
{
    Vector2 k1v1 = ComputeAcceleration(body1, body2);
    Vector2 k1p1 = body1->velocity;

    Vector2 k1v2 = ComputeAcceleration(body2, body1);
    Vector2 k1p2 = body2->velocity;

    PointMass mid1 = {
        { body1->position.x + 0.5f * k1p1.x * dt, body1->position.y + 0.5f * k1p1.y * dt },
        { body1->velocity.x + 0.5f * k1v1.x * dt, body1->velocity.y + 0.5f * k1v1.y * dt },
        body1->mass
    };

    PointMass mid2 = {
        { body2->position.x + 0.5f * k1p2.x * dt, body2->position.y + 0.5f * k1p2.y * dt },
        { body2->velocity.x + 0.5f * k1v2.x * dt, body2->velocity.y + 0.5f * k1v2.y * dt },
        body2->mass
    };

    Vector2 k2v1 = ComputeAcceleration(&mid1, &mid2);
    Vector2 k2p1 = mid1.velocity;

    Vector2 k2v2 = ComputeAcceleration(&mid2, &mid1);
    Vector2 k2p2 = mid2.velocity;

    body1->velocity.x += dt * k2v1.x;
    body1->velocity.y += dt * k2v1.y;
    body1->position.x += dt * k2p1.x;
    body1->position.y += dt * k2p1.y;

    body2->velocity.x += dt * k2v2.x;
    body2->velocity.y += dt * k2v2.y;
    body2->position.x += dt * k2p2.x;
    body2->position.y += dt * k2p2.y;
}

void UpdateRK3_Classic(PointMass* body1, PointMass* body2, const float dt)
{
    Vector2 k1v1 = ComputeAcceleration(body1, body2);
    Vector2 k1p1 = body1->velocity;

    Vector2 k1v2 = ComputeAcceleration(body2, body1);
    Vector2 k1p2 = body2->velocity;

    PointMass mid1 = {
        { body1->position.x + 0.5f * k1p1.x * dt, body1->position.y + 0.5f * k1p1.y * dt },
        { body1->velocity.x + 0.5f * k1v1.x * dt, body1->velocity.y + 0.5f * k1v1.y * dt },
        body1->mass
    };

    PointMass mid2 = {
        { body2->position.x + 0.5f * k1p2.x * dt, body2->position.y + 0.5f * k1p2.y * dt },
        { body2->velocity.x + 0.5f * k1v2.x * dt, body2->velocity.y + 0.5f * k1v2.y * dt },
        body2->mass
    };

    Vector2 k2v1 = ComputeAcceleration(&mid1, &mid2);
    Vector2 k2p1 = mid1.velocity;

    Vector2 k2v2 = ComputeAcceleration(&mid2, &mid1);
    Vector2 k2p2 = mid2.velocity;

    PointMass end1 = {
        { body1->position.x + dt * (-k1p1.x + 2.0f * k2p1.x),
          body1->position.y + dt * (-k1p1.y + 2.0f * k2p1.y) },
        { body1->velocity.x + dt * (-k1v1.x + 2.0f * k2v1.x),
          body1->velocity.y + dt * (-k1v1.y + 2.0f * k2v1.y) },
        body1->mass
    };

    PointMass end2 = {
        { body2->position.x + dt * (-k1p2.x + 2.0f * k2p2.x),
          body2->position.y + dt * (-k1p2.y + 2.0f * k2p2.y) },
        { body2->velocity.x + dt * (-k1v2.x + 2.0f * k2v2.x),
          body2->velocity.y + dt * (-k1v2.y + 2.0f * k2v2.y) },
        body2->mass
    };

    Vector2 k3v1 = ComputeAcceleration(&end1, &end2);
    Vector2 k3p1 = end1.velocity;

    Vector2 k3v2 = ComputeAcceleration(&end2, &end1);
    Vector2 k3p2 = end2.velocity;

    body1->velocity.x += (dt / 6.0f) * (k1v1.x + 4.0f * k2v1.x + k3v1.x);
    body1->velocity.y += (dt / 6.0f) * (k1v1.y + 4.0f * k2v1.y + k3v1.y);
    body1->position.x += (dt / 6.0f) * (k1p1.x + 4.0f * k2p1.x + k3p1.x);
    body1->position.y += (dt / 6.0f) * (k1p1.y + 4.0f * k2p1.y + k3p1.y);

    body2->velocity.x += (dt / 6.0f) * (k1v2.x + 4.0f * k2v2.x + k3v2.x);
    body2->velocity.y += (dt / 6.0f) * (k1v2.y + 4.0f * k2v2.y + k3v2.y);
    body2->position.x += (dt / 6.0f) * (k1p2.x + 4.0f * k2p2.x + k3p2.x);
    body2->position.y += (dt / 6.0f) * (k1p2.y + 4.0f * k2p2.y + k3p2.y);
}

void UpdateRK4(PointMass* body1, PointMass* body2, const float dt)
{
    Vector2 k1v1 = ComputeAcceleration(body1, body2);
    Vector2 k1p1 = body1->velocity;
    Vector2 k1v2 = ComputeAcceleration(body2, body1);
    Vector2 k1p2 = body2->velocity;

    Vector2 midVelocity1 = { body1->velocity.x + 0.5f * k1v1.x * dt, body1->velocity.y + 0.5f * k1v1.y * dt };
    Vector2 midPosition1 = { body1->position.x + 0.5f * k1p1.x * dt, body1->position.y + 0.5f * k1p1.y * dt };
    Vector2 midVelocity2 = { body2->velocity.x + 0.5f * k1v2.x * dt, body2->velocity.y + 0.5f * k1v2.y * dt };
    Vector2 midPosition2 = { body2->position.x + 0.5f * k1p2.x * dt, body2->position.y + 0.5f * k1p2.y * dt };

    PointMass mass1 = (PointMass){ midPosition1, midVelocity1, body1->mass };
    PointMass mass2 = (PointMass){ midPosition2, midVelocity2, body2->mass };

    Vector2 k2v1 = ComputeAcceleration(&mass1, &mass2);
    Vector2 k2p1 = midVelocity1;
    Vector2 k2v2 = ComputeAcceleration(&mass2, &mass1);
    Vector2 k2p2 = midVelocity2;

    Vector2 midVelocity1_2 = { body1->velocity.x + 0.5f * k2v1.x * dt, body1->velocity.y + 0.5f * k2v1.y * dt };
    Vector2 midPosition1_2 = { body1->position.x + 0.5f * k2p1.x * dt, body1->position.y + 0.5f * k2p1.y * dt };
    Vector2 midVelocity2_2 = { body2->velocity.x + 0.5f * k2v2.x * dt, body2->velocity.y + 0.5f * k2v2.y * dt };
    Vector2 midPosition2_2 = { body2->position.x + 0.5f * k2p2.x * dt, body2->position.y + 0.5f * k2p2.y * dt };

    mass1 = (PointMass){ midPosition1_2, midVelocity1_2, body1->mass };
    mass2 = (PointMass){ midPosition2_2, midVelocity2_2, body2->mass };

    Vector2 k3v1 = ComputeAcceleration(&mass1, &mass2);
    Vector2 k3p1 = midVelocity1_2;
    Vector2 k3v2 = ComputeAcceleration(&mass2, &mass1);
    Vector2 k3p2 = midVelocity2_2;

    Vector2 endVelocity1 = { body1->velocity.x + k3v1.x * dt, body1->velocity.y + k3v1.y * dt };
    Vector2 endPosition1 = { body1->position.x + k3p1.x * dt, body1->position.y + k3p1.y * dt };
    Vector2 endVelocity2 = { body2->velocity.x + k3v2.x * dt, body2->velocity.y + k3v2.y * dt };
    Vector2 endPosition2 = { body2->position.x + k3p2.x * dt, body2->position.y + k3p2.y * dt };

    mass1 = (PointMass){ endPosition1, endVelocity1, body1->mass };
    mass2 = (PointMass){ endPosition2, endVelocity2, body2->mass };

    Vector2 k4v1 = ComputeAcceleration(&mass1, &mass2);
    Vector2 k4p1 = endVelocity1;
    Vector2 k4v2 = ComputeAcceleration(&mass2, &mass1);
    Vector2 k4p2 = endVelocity2;

    body1->velocity.x += (dt / 6.0f) * (k1v1.x + 2 * k2v1.x + 2 * k3v1.x + k4v1.x);
    body1->velocity.y += (dt / 6.0f) * (k1v1.y + 2 * k2v1.y + 2 * k3v1.y + k4v1.y);
    body1->position.x += (dt / 6.0f) * (k1p1.x + 2 * k2p1.x + 2 * k3p1.x + k4p1.x);
    body1->position.y += (dt / 6.0f) * (k1p1.y + 2 * k2p1.y + 2 * k3p1.y + k4p1.y);

    body2->velocity.x += (dt / 6.0f) * (k1v2.x + 2 * k2v2.x + 2 * k3v2.x + k4v2.x);
    body2->velocity.y += (dt / 6.0f) * (k1v2.y + 2 * k2v2.y + 2 * k3v2.y + k4v2.y);
    body2->position.x += (dt / 6.0f) * (k1p2.x + 2 * k2p2.x + 2 * k3p2.x + k4p2.x);
    body2->position.y += (dt / 6.0f) * (k1p2.y + 2 * k2p2.y + 2 * k3p2.y + k4p2.y);
}

float ComputeKineticEnergy(const PointMass* body)
{
    const float speedSquared = body->velocity.x * body->velocity.x + body->velocity.y * body->velocity.y;
    return 0.5f * body->mass * speedSquared;
}

float ComputePotentialEnergy(const PointMass* a, const PointMass* b)
{
    const Vector2 direction = { b->position.x - a->position.x, b->position.y - a->position.y };
    float distance = sqrtf(direction.x * direction.x + direction.y * direction.y);
    if (distance < 1.0f) distance = 1.0f; // Avoid division by zero
    return -(G * a->mass * b->mass) / distance;
}

void FormatFloat(char* buffer, const size_t size, const float value, const int precision)
{
    snprintf(buffer, size, "%.*f", precision, value);
}

char kineticStr[64];
char potentialStr[64];
char totalStr[64];
char text1[128];
char text2[128];
char text3[128];

int main()
{
    const int screenWidth = 800;
    const int screenHeight = 600;

    InitWindow(screenWidth, screenHeight, "2D Physics Examples");
    SetTargetFPS(60);

    const Vector2 centerMass = (Vector2){400, 300};

    PointMass body1 = (PointMass){ { centerMass.x - 100, centerMass.y }, { 0, 0 }, 10.0f };
    PointMass body2 = (PointMass){ { centerMass.x + 100, centerMass.y }, { 0, 0 }, 10.0f };

    const float distance1 = sqrtf(powf(body2.position.x - body1.position.x, 2) + powf(body2.position.y - body1.position.y, 2));
    const float orbitalSpeed1 = sqrtf(G * (body1.mass + body2.mass) / distance1);

    body1.velocity = (Vector2){ 0, -orbitalSpeed1 * (body2.mass / (body1.mass + body2.mass)) };
    body2.velocity = (Vector2){ 0, orbitalSpeed1 * (body1.mass / (body1.mass + body2.mass)) };

    while (!WindowShouldClose()) 
    {
        if (IsKeyPressed(KEY_ONE)) currentMethod = RK1;
        if (IsKeyPressed(KEY_TWO)) currentMethod = RK2;
        if (IsKeyPressed(KEY_THREE)) currentMethod = RK3;
        if (IsKeyPressed(KEY_FOUR)) currentMethod = RK4;

        if (currentMethod == RK1) UpdateRK1(&body1, &body2, TIME_STEP);
        else if (currentMethod == RK2) UpdateRK2_Midpoint(&body1, &body2, TIME_STEP);
        else if (currentMethod == RK3) UpdateRK3_Classic(&body1, &body2, TIME_STEP);
        else if (currentMethod == RK4) UpdateRK4(&body1, &body2, TIME_STEP);

        const float kineticEnergy = ComputeKineticEnergy(&body1) + ComputeKineticEnergy(&body2);
        const float potentialEnergy = ComputePotentialEnergy(&body1, &body2);
        const float totalEnergy = kineticEnergy + potentialEnergy;

        BeginDrawing();
        ClearBackground(RAYWHITE);

        DrawCircleV(body1.position, 10, RED);
        DrawCircleV(body2.position, 10, BLUE);

        FormatFloat(kineticStr, sizeof(kineticStr), kineticEnergy, 3);
        FormatFloat(potentialStr, sizeof(potentialStr), potentialEnergy, 3);
        FormatFloat(totalStr, sizeof(totalStr), totalEnergy, 3);

        snprintf(text1, sizeof(text1), "Kinetic Energy: %s", kineticStr);
        snprintf(text2, sizeof(text2), "Potential Energy: %s", potentialStr);
        snprintf(text3, sizeof(text3), "Total Energy: %s", totalStr);

        DrawText(currentMethod == RK1 ? "Method: Euler" : "Method: RK4", 10, 10, 20, BLACK);
        DrawText("Press 1 for Euler, 2 for RK4", 10, 40, 20, BLACK);
        DrawText(text1, 10, 80, 20, BLACK);
        DrawText(text2, 10, 110, 20, BLACK);
        DrawText(text3, 10, 140, 20, BLACK);

        DrawText("BACKSPACE to return, R to reset positions", 10, GetScreenHeight() - 25, 20, BLACK);
        EndDrawing();
    }

    CloseWindow();
}