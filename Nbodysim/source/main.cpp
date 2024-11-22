/************************************************************
*  Author:         Nicholas Tibbetts
*  Date:           
*  Course Code:    
*  License:        Copyright 2024 Nic Tibbetts
*  References:     _
*  Worked with:    _
*  Description:    _
***********************************************************/

#include "Simulation.hpp"
#include "raylib.h"
#include <thread>
#include <atomic>
#include <mutex>
#include <memory>

std::atomic<bool> PAUSED{false};
std::mutex UPDATE_LOCK;
std::vector<Body> SHARED_BODIES;
std::vector<Node> SHARED_QUADTREE;
std::vector<Body> SPAWN_QUEUE;


struct SystemMetrics
{
    float centralMass;
    float totalMass;
    float totalKineticEnergy;
    float totalPotentialEnergy;
    float averageOrbitalPeriod;
    float netForce;
    float averageSpeed;
};


SystemMetrics calculateMetrics(const std::vector<Body>& bodies) {
    SystemMetrics metrics = {};
    
    // "finds" central mass (ie body closest to center of mass)
    Vec2 centerOfMass = Vec2::zero();
    float totalMass = 0;
    
    // calculates the center of mass first
    for (const auto& body : bodies) {
        centerOfMass += body.pos * body.mass;
        totalMass += body.mass;
    }
    centerOfMass = centerOfMass / totalMass;
    
    // then finds body closest to center of mass with significant mass
    const Body* centralBody = nullptr;
    float minDist = std::numeric_limits<float>::max();
    const float MASS_THRESHOLD = totalMass * 0.1f; // 10% of total mass
    
    for (const auto& body : bodies) {
        if (body.mass > MASS_THRESHOLD) {
            float dist = (body.pos - centerOfMass).mag_sq();
            if (dist < minDist) {
                minDist = dist;
                centralBody = &body;
            }
        }
    }
    
    if (!centralBody) return metrics;
    
    metrics.centralMass = centralBody->mass;
    metrics.totalMass = totalMass;
    
    // will count bodies in stable orbits
    size_t stableBodyCount = 0;
    // NOTE: these need to be dynamic change based on your scale instead of manually set
    const float ESCAPE_VELOCITY_SQ = 2.0f; // adjust this threshold
    const float MAX_ORBITAL_RADIUS = 1000.0f; // adjust based on the scale
    
    for (const auto& body : bodies) {
        if (&body != centralBody) {
            Vec2 r = body.pos - centralBody->pos;
            float dist = r.mag();
            
            // skips ejected/unstable bodies
            if (dist > MAX_ORBITAL_RADIUS) continue;
            
            Vec2 relativeVel = body.vel - centralBody->vel;
            float speedSq = relativeVel.mag_sq();
            float escapeSpeedSq = 2.0f * centralBody->mass / dist;
            
            // it will only include bodies in stable orbits
            if (speedSq < escapeSpeedSq * ESCAPE_VELOCITY_SQ && dist < MAX_ORBITAL_RADIUS) {
                stableBodyCount++;
                
                // kinetic Energy relative to central body
                metrics.totalKineticEnergy += 0.5f * body.mass * speedSq;
                
                // potential Energy with corrected sign
                metrics.totalPotentialEnergy += -1.0f * body.mass * centralBody->mass / dist;
                
                // orbital period using corrected velocity
                float semiMajorAxis = dist; // simplified NOTE: should use actual orbital elements
                float period = 2.0f * PI * std::sqrt(semiMajorAxis * semiMajorAxis * semiMajorAxis / 
                                                   (centralBody->mass));
                metrics.averageOrbitalPeriod += period;
                
                // Force(net)  (gravitational)
                metrics.netForce += body.mass * centralBody->mass / (dist * dist);
                
                // average speed relative to central body
                metrics.averageSpeed += std::sqrt(speedSq);
            }
        }
    }
    
    // lastly calculate the averages only for stable bodies
    if (stableBodyCount > 0) {
        metrics.averageOrbitalPeriod /= stableBodyCount;
        metrics.averageSpeed /= stableBodyCount;
    }
    
    return metrics;
}

// converts simulation coordinates to screen coordinates
Vector2 worldToScreen(Vec2 worldPos, float scale, Vector2 center)
{
    return Vector2{
        (worldPos.x * scale) + GetScreenWidth() / 20.0f + center.x,
        (worldPos.y * scale) + GetScreenHeight() / 20.0f + center.y};
}

void drawQuadtreeNode(const Node &node, const std::vector<Node> &nodes, float scale, Vector2 center)
{
    Vector2 pos = worldToScreen(node.pos, scale, center);
    DrawCircleV(pos, 2.0f, RED);
    float size = node.quad.size * scale;
    Vector2 quadPos = worldToScreen(node.quad.center, scale, center);

    // node boundary
    DrawRectangleLinesEx(
        {quadPos.x - size / 2, quadPos.y - size / 2, size, size},
        1,
        {100, 100, 100, 100});

    // children if this is a branch node
    if (node.is_branch())
    {
        for (size_t i = 0; i < 4; ++i)
        {
            drawQuadtreeNode(nodes[node.children + i], nodes, scale, center);
        }
    }
}
void drawBlackHole(const Body &body, float scale, Vector2 center)
{
    Vector2 screenPos = worldToScreen(body.pos, scale, center);
    float screenRadius = std::max(2.0f, body.radius * scale);

    // Gravitational lensing effect
    for (int i = 4; i >= 0; i--)
    {
        float alpha = (1.0f - i / 4.0f) * 1.1f;
        DrawCircleGradient(
            screenPos.x, screenPos.y,
            screenRadius * (1.0f + i * 1.4f),
            (Color){255, 255, 237, (unsigned char)(alpha * 255)},
            (Color){0, 0, 0, 0});
    }

    // accretion disk
    float diskStartRadius = screenRadius * 2.1f;
    float diskEndRadius = screenRadius * 10.51f;
    int segments = 10096;

    for (int i = 0; i < segments; i++)
    {
        float angle = (i * 300.0f / segments) * DEG2RAD;
        float nextAngle = ((i + 1) * 390.0f / segments) * DEG2RAD;

        // an attempt at a relativistic beaming - brighter on approaching side
        float brightness = 1.4f + 1.0f * (10.5f + cosf(angle));
        // a vertical distortion from gravitational lensing
        float distortion = 0.55f + 0.10f * (1.02f - tanf(angle * 12.0f));

        Vector2 p1 = {
            screenPos.x + cosf(angle) * diskStartRadius,
            screenPos.y + sinf(angle) * diskStartRadius * distortion};
        Vector2 p2 = {
            screenPos.x + cosf(nextAngle) * diskStartRadius,
            screenPos.y + sinf(nextAngle) * diskStartRadius * distortion};
        Vector2 p3 = {
            screenPos.x + cosf(angle) * diskEndRadius,
            screenPos.y + sinf(angle) * diskEndRadius * distortion};
        Vector2 p4 = {
            screenPos.x + cosf(nextAngle) * diskEndRadius,
            screenPos.y + sinf(nextAngle) * diskEndRadius * distortion};

        // disk segments with color gradient
        Color diskColor = {
            (unsigned char)(3 * brightness), // R
            (unsigned char)(2 * brightness), // G
            (unsigned char)(6 * brightness), // B
            (unsigned char)(2)               // higher alpha for visibility
        };

        DrawTriangle(p1, p3, p2, diskColor);
        DrawTriangle(p2, p3, p4, diskColor);
    }

    // event horizon (pure black with subtle blue edge)
    DrawCircleGradient(
        screenPos.x, screenPos.y,
        screenRadius * 1.03f,
        (Color){0, 0, 40, 255},
        BLACK);
    DrawCircleV(screenPos, screenRadius, BLACK);

    // einstein ring (photon sphere)
    float ringRadius = screenRadius;
    float ringThickness = screenRadius * 0.011f;
    DrawRing(
        screenPos,
        ringRadius - ringThickness / 2,
        ringRadius + ringThickness / 2,
        0, 360,
        60, // segments
        (Color){255, 225, 210, 255});
}

inline Color getStarColorWithBrightness(const Body &body, float brightness = 1.0f)
{
    Color base;
    float mass = body.mass;

    // color ranges and thresholds (completely arbitrary)
    const Color colors[] = {
        {200, 0, 0, 255},     // deep red (Brown Dwarfs)
        {255, 50, 0, 255},    // orange red (Red Dwarfs)
        {255, 100, 0, 255},   // deep orange (Orange Dwarfs)
        {255, 150, 50, 255},  // light orange (Transition to Yellow)
        {255, 240, 150, 255}, // yellow (Sun-like stars)
        {255, 255, 200, 255}, // light yellow (Transition to White)
        {219, 233, 244, 255}, // bluish white (White Stars)
        {173, 216, 230, 255}, // light blue (Blue-White Stars)
        {100, 100, 255, 255}, // blue (Blue Stars)
        {0, 0, 255, 255},     // deep blue (Hyper-giant Blue Stars)
        {5, 5, 5, 2}          // invisible/black (neutron stars or potentially extreme cases)
    };

    const float thresholds[] = {0.08f, 0.4f, 0.8f, 1.2f, 1.5f, 2.5f, 5.0f, 15.0f, 25.0f, 50.0f};

    // essentially interpolates color between two thresholds if the mass is between them
    // base = colors[sizeof(thresholds) / sizeof(float)]; // defaults to last color
    // for (size_t i = 0; i < sizeof(thresholds) / sizeof(float) - 1; i++) {
    //     if (mass < thresholds[i + 1]) {
    //         float t = (mass - thresholds[i]) / (thresholds[i + 1] - thresholds[i]);
    //         base.r = (unsigned char)(colors[i].r * (1 - t) + colors[i + 1].r * t);
    //         base.g = (unsigned char)(colors[i].g * (1 - t) + colors[i + 1].g * t);
    //         base.b = (unsigned char)(colors[i].b * (1 - t) + colors[i + 1].b * t);
    //         break;
    //     }
    // }

    base = colors[sizeof(thresholds) / sizeof(float)]; // defaults to last color
    for (size_t i = 0; i < sizeof(thresholds) / sizeof(float); i++)
    {
        if (mass < thresholds[i])
        {
            base = colors[i];
            break;
        }
    }

    // brightness adjustment
    return (Color){
        (unsigned char)(base.r * brightness),
        (unsigned char)(base.g * brightness),
        (unsigned char)(base.b * brightness),
        255};
}

void simulation_thread(std::shared_ptr<Simulation> simulation)
{
    while (!WindowShouldClose())
    {
        if (PAUSED)
        {
            std::this_thread::yield();
            continue;
        }
        simulation->step();

        {
            std::lock_guard<std::mutex> lock(UPDATE_LOCK);
            // update the shared state
            SHARED_BODIES = simulation->bodies;
            SHARED_QUADTREE = simulation->quadtree.nodes;
        }

        // sleep for a short duration to prevent maxing out CPU
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
}

int main()
{
    const int WIDTH = 1200;
    const int HEIGHT = 900;
    // float scale = 0.5f;
    // Vector2 center = {WIDTH / 2.0f, HEIGHT / 2.0f};
    float scale = 1.0f;
    Vector2 center = {0, 0};
    bool showQuadtree = false;

    InitWindow(WIDTH, HEIGHT, "N-Body Simulation");
    // SetTargetFPS(60);
    // SetTargetFPS(GetFPS() < 30 ? 60 : 0);  // 0 means uncapped

    auto simulation = std::make_shared<Simulation>();

    // starts simulation in separate thread
    std::thread sim_thread(simulation_thread, simulation);
    sim_thread.detach();

    Camera2D camera = {
        .offset = {0, 0},
        .target = {0, 0},
        .rotation = 0,
        .zoom = 1.0f};
    camera.target = {0, 0};
    camera.offset = {WIDTH / 2.0f, HEIGHT / 2.0f};
    camera.rotation = 0.0f;
    camera.zoom = 1.0f;

    while (!WindowShouldClose())
    {
        if (IsKeyPressed(KEY_SPACE))
        {
            PAUSED = !PAUSED;
        }
        if (IsKeyPressed(KEY_Q))
        {
            showQuadtree = !showQuadtree;
        }

        if (IsKeyDown(KEY_W))
            camera.offset.y += 10.0f;
        if (IsKeyDown(KEY_S))
            camera.offset.y -= 10.0f;
        if (IsKeyDown(KEY_A))
            camera.offset.x += 10.0f;
        if (IsKeyDown(KEY_D))
            camera.offset.x -= 10.0f;
        if (IsKeyDown(KEY_R))
            scale *= 1.02f;
        if (IsKeyDown(KEY_F))
            scale *= 0.98f;

        BeginDrawing();
        ClearBackground(BLACK);

        BeginMode2D(camera);

        // draw bodies
        {
            std::lock_guard<std::mutex> lock(UPDATE_LOCK);

            // Draw quadtree if enabled
            if (showQuadtree && !SHARED_QUADTREE.empty())
            {
                drawQuadtreeNode(SHARED_QUADTREE[0], SHARED_QUADTREE, scale, center);
            }

            // finds the largest body (assumed to be the central mass)
            const Body *centralMass = nullptr;
            float maxRadius = 0;

            for (const auto &body : SHARED_BODIES)
            {
                if (body.radius > maxRadius)
                {
                    maxRadius = body.radius;
                    centralMass = &body;
                }
            }

            // regular bodies
            for (const auto &body : SHARED_BODIES)
            {
                if (&body == centralMass)
                    continue;

                Vector2 screenPos = worldToScreen(body.pos, scale, center);
                float screenRadius = std::max(1.0f, body.radius * scale);

                // Skip tiny stars for performance
                if (screenRadius < 1.0f)
                    continue;

                // Single draw call with proper color
                DrawCircleV(screenPos, screenRadius, getStarColorWithBrightness(body, 3.0f));
            }

            // draw central mass
            if (centralMass)
            {
                drawBlackHole(*centralMass, scale, center);
            }
        }

        EndMode2D();

        DrawFPS(10, 100);

        // UI goodies
        DrawText("Bodies: ", 10, 10, 20, WHITE);
        DrawText(TextFormat(" %zu", SHARED_BODIES.size()), 80, 10, 20, GetFPS() >= 30 ? GREEN : (GetFPS() >= 15 ? ORANGE : RED));
        DrawText(TextFormat("Scale: %.2f", (double)scale), 10, 35, 20, WHITE);
        DrawText(PAUSED ? "PAUSED" : "RUNNING", 10, 60, 20, PAUSED ? RED : GREEN);
        DrawText("Controls:", 10, HEIGHT - 100, 20, WHITE);
        DrawText("WASD: Pan | R/F: Zoom | Space: Pause | Q: Toggle Quadtree", 10, HEIGHT - 75, 20, WHITE);

        {
            std::lock_guard<std::mutex> lock(UPDATE_LOCK);
            SystemMetrics metrics = calculateMetrics(SHARED_BODIES);

            // Displays all of the metrics in top right corner
            int rightX = GetScreenWidth() - 150;
            int y = 10;
            int lineHeight = 20;

            DrawText("System Metrics:", rightX, y, 15, WHITE);
            y += lineHeight * 1.5f;

            DrawText(TextFormat("Central Mass: %.2e", (double)metrics.centralMass), rightX, y, 5, WHITE);
            y += lineHeight;

            DrawText(TextFormat("Total Mass: %.2e", (double)metrics.totalMass), rightX, y, 5, WHITE);
            y += lineHeight;

            DrawText(TextFormat("Avg Speed: %.2f", (double)metrics.averageSpeed), rightX, y, 5, WHITE);
            y += lineHeight;

            DrawText(TextFormat("Kinetic Energy: %.2e", (double)metrics.totalKineticEnergy), rightX, y, 5, WHITE);
            y += lineHeight;

            DrawText(TextFormat("Potential Energy: %.2e", (double)metrics.totalPotentialEnergy), rightX, y, 5, WHITE);
            y += lineHeight;

            float totalEnergy = metrics.totalKineticEnergy + metrics.totalPotentialEnergy;
            DrawText(TextFormat("Total Energy: %.2e", (double)totalEnergy), rightX, y, 5, WHITE);
            y += lineHeight;

            DrawText(TextFormat("Avg Orbital Period: %.2f", (double)metrics.averageOrbitalPeriod), rightX, y, 5, WHITE);
            y += lineHeight;

            DrawText(TextFormat("Net Force: %.2e", (double)metrics.netForce), rightX, y, 5, WHITE);
        }

        EndDrawing();
    }

    CloseWindow();
    return 0;
}