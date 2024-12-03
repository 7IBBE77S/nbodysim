/************************************************************
 *  Author:         Nicholas Tibbetts
 *  Date:
 *  Course Code:
 *  License:        Copyright 2024 Nic Tibbetts
 *  References:     _
 *  Worked with:    _
 *  Description:    _
 ***********************************************************/

#include "raylib.h"
#include "rlgl.h"
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wall"
#pragma clang diagnostic ignored "-Weverything"

#define RAYGUI_IMPLEMENTATION
#include "raygui.h"

#pragma clang diagnostic pop

#include "Simulation.hpp"
#include "Vec2.hpp"

#include <thread>
#include <mutex>
#include <memory>
#include <unordered_map>
#include <thread>
#include <future>

std::atomic<bool> PAUSED{false};
std::atomic<bool> PERFORMANCE_MODE{false};
std::atomic<bool> SHOW_CONNECTIONS{false};
std::atomic<bool> SHOW_QUADTREE{false};
std::atomic<bool> SHOW_BODIES{true};

std::mutex UPDATE_LOCK;
std::atomic<float> SIMULATION_DT{0.01f};
std::vector<Body> SHARED_BODIES;
std::vector<Node> SHARED_QUADTREE;
std::vector<Body> SPAWN_QUEUE;

std::condition_variable frameCondition;
std::mutex simulationMutex;
bool frameRendered = false;

constexpr int MAX_BODIES_PER_CELL = 16;

constexpr float MAX_DISTANCE = 1000.0f;
constexpr int MAX_CONNECTIONS = 5;

constexpr float SLIDER_MIN_DT = 0.001f;
constexpr float SLIDER_MAX_DT = 0.1f;

struct GridCell
{
    const Body *bodies[MAX_BODIES_PER_CELL];
    int count = 0;
};

bool firstFrame = true;
std::unordered_map<uint64_t, GridCell> spatialGrid;

constexpr float MAX_GRID_DISTANCE = MAX_DISTANCE * 100.0f;

uint64_t getGridKey(const Vec2 &pos)
{
    // clamp positions to prevent grid explosion
    float x = std::clamp(pos.x, -MAX_GRID_DISTANCE, MAX_GRID_DISTANCE);
    float y = std::clamp(pos.y, -MAX_GRID_DISTANCE, MAX_GRID_DISTANCE);

    constexpr int64_t GRID_OFFSET = 1000000;
    int64_t gridX = static_cast<int64_t>(x / MAX_DISTANCE) + GRID_OFFSET;
    int64_t gridY = static_cast<int64_t>(y / MAX_DISTANCE) + GRID_OFFSET;

    return (static_cast<uint64_t>(gridX) << 32) | static_cast<uint64_t>(gridY);
}

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

SystemMetrics calculateMetrics(const std::vector<Body> &bodies)
{
    SystemMetrics metrics = {};
    static float smoothedDt = SIMULATION_DT.load();
    smoothedDt = smoothedDt * 0.95f + SIMULATION_DT.load() * 0.05f;

    Vec2 minBound = Vec2{std::numeric_limits<float>::max(), std::numeric_limits<float>::max()};
    Vec2 maxBound = Vec2{std::numeric_limits<float>::lowest(), std::numeric_limits<float>::lowest()};

    for (const auto &body : bodies)
    {
        minBound.x = std::min(minBound.x, body.pos.x);
        minBound.y = std::min(minBound.y, body.pos.y);
        maxBound.x = std::max(maxBound.x, body.pos.x);
        maxBound.y = std::max(maxBound.y, body.pos.y);
    }

    Vec2 systemSize = maxBound - minBound;
    float systemRadius = std::sqrt(systemSize.mag_sq()) * 0.5f;
    const float BASE_ORBITAL_RADIUS = std::max(1000.0f, systemRadius);

    Vec2 centerOfMass = Vec2::zero();
    float totalMass = 0;

    for (const auto &body : bodies)
    {
        centerOfMass += body.pos * body.mass;
    }
    centerOfMass = centerOfMass / totalMass;

    const Body *centralBody = nullptr;
    float minDist = std::numeric_limits<float>::max();
    const float MASS_THRESHOLD = totalMass * 0.1f;

    for (const auto &body : bodies)
    {
        if (body.mass > MASS_THRESHOLD)
        {
            float dist = (body.pos - centerOfMass).mag_sq();
            if (dist < minDist)
            {
                minDist = dist;
                centralBody = &body;
            }
        }
    }

    if (!centralBody)
        return metrics;

    metrics.centralMass = centralBody->mass;
    metrics.totalMass = totalMass;

    size_t stableBodyCount = 0;
    const float BASE_ESCAPE_VELOCITY_SQ = 2.0f;

    float escapeThreshold = BASE_ESCAPE_VELOCITY_SQ * (1.0f + std::log10(smoothedDt + 1.0f));

    for (const auto &body : bodies)
    {
        if (&body != centralBody)
        {
            Vec2 r = body.pos - centralBody->pos;
            float dist = r.mag();

            if (dist > BASE_ORBITAL_RADIUS * 2.0f)
                continue;

            Vec2 relativeVel = body.vel - centralBody->vel;
            float speedSq = relativeVel.mag_sq();

            float escapeSpeedSq = 2.0f * centralBody->mass / dist;

            if (speedSq < escapeSpeedSq * escapeThreshold)
            {
                stableBodyCount++;

                metrics.totalKineticEnergy += 0.5f * body.mass * speedSq;
                metrics.totalPotentialEnergy += -1.0f * body.mass * centralBody->mass / dist;

                float semiMajorAxis = dist;
                float period = 2.0f * PI * std::sqrt(semiMajorAxis * semiMajorAxis * semiMajorAxis / centralBody->mass);
                metrics.averageOrbitalPeriod += period;

                metrics.netForce += body.mass * centralBody->mass / (dist * dist);
                metrics.averageSpeed += std::sqrt(speedSq);
            }
        }
    }

    if (stableBodyCount > 0)
    {
        metrics.averageOrbitalPeriod /= stableBodyCount;
        metrics.averageSpeed /= stableBodyCount;
    }

    metrics.totalKineticEnergy *= smoothedDt;
    metrics.totalPotentialEnergy *= smoothedDt;
    metrics.netForce *= smoothedDt;
    metrics.averageSpeed *= std::sqrt(smoothedDt);
    metrics.averageOrbitalPeriod *= smoothedDt;

    return metrics;
}

Vector2 worldToScreen(Vec2 worldPos, float scale, Vector2 center)
{
    return Vector2{
        (worldPos.x - center.x) * scale + GetScreenWidth() / 2.0f,
        (worldPos.y - center.y) * scale + GetScreenHeight() / 2.0f};
}

bool isInScreenBounds(const Vec2 &worldPos, float scale, Vector2 center, float margin = 1000.0f)
{
    float scaledMargin = margin / scale;

    Vector2 screenPos = worldToScreen(worldPos, scale, center);
    float width = GetScreenWidth() + scaledMargin;
    float height = GetScreenHeight() + scaledMargin;

    return screenPos.x >= -scaledMargin && screenPos.x <= width &&
           screenPos.y >= -scaledMargin && screenPos.y <= height;
}

struct RenderStats
{
    float lastFrameTime = 0.0f;
    float targetFrameTime = 1.0f / 60.0f; // target  FPS
    // float targetFrameTime = 0.0f; // Uncapped FPS
    int totalVisibleBodies = 0;
    int totalConnections = 0;
    float smoothedVisibleRatio = 1.0f;
    float smoothedPerformanceRatio = 1.0f;
} renderStats;

struct Connection
{
    Vector2 pos1;
    Vector2 pos2;
    unsigned char alpha;
};
// granular connection mode to the clustered cell mode
void drawConnections(const std::vector<Body> &bodies, float scale, Vector2 center)
{
    if (!SHOW_CONNECTIONS)
        return;

    float frameStart = GetTime();
    spatialGrid.clear();

    // inverse relationship - fewer connections when zoomed out*
    float zoomFactor = std::max(0.1f, scale);
    int gridLevel = std::max(0, static_cast<int>(-std::log2(zoomFactor))); // more levels when zoomed out
    float baseDistance = MAX_DISTANCE;
    float adaptiveDistance = baseDistance * (1.0f / zoomFactor); // increases with zoom out
    float adaptiveDistanceSq = adaptiveDistance * adaptiveDistance;

    // Connection count decreases with zoom out
    int adaptiveMaxConnections = static_cast<int>(MAX_CONNECTIONS * (1.0f / zoomFactor));

    // Keep alpha more visible
    float baseAlpha = 255.0f;
    float adaptiveAlpha = std::max(50.0f, baseAlpha * zoomFactor); // Minimum alpha of 50

    // Populate spatial grid without culling
    for (const auto &body : bodies)
    {
        uint64_t key = getGridKey(body.pos);
        auto &cell = spatialGrid[key];
        if (cell.count < MAX_BODIES_PER_CELL)
        {
            cell.bodies[cell.count++] = &body;
        }
    }

    // Batch render connections
    rlBegin(RL_LINES);
    int totalConnections = 0;

    for (const auto &[key, cell] : spatialGrid)
    {
        if (cell.count == 0)
            continue;
        if (gridLevel > 2 && cell.count > 1)
        {
            //     // Draw one connection per cell to represent the group
            //     Vec2 cellCenter = Vec2::zero();
            //     for (int i = 0; i < cell.count; i++)
            //     {
            //         cellCenter = cellCenter + cell.bodies[i]->pos;
            //     }
            //     cellCenter = cellCenter * (1.0f / cell.count);

            //     Vector2 screenPos = worldToScreen(cellCenter, scale, center);

            //     // Draw connections to neighboring cells only
            //     for (int dx = -1; dx <= 1; dx++)
            //     {
            //         for (int dy = -1; dy <= 1; dy++)
            //         {
            //             if (dx == 0 && dy == 0)
            //                 continue;

            //             uint64_t nKey = getGridKey(Vec2{
            //                 cellCenter.x + dx * baseDistance,
            //                 cellCenter.y + dy * baseDistance});

            //             const auto &nCellIt = spatialGrid.find(nKey);
            //             if (nCellIt != spatialGrid.end() && nCellIt->second.count > 0)
            //             {
            //                 Vec2 nCellCenter = Vec2::zero();
            //                 for (int i = 0; i < nCellIt->second.count; i++)
            //                 {
            //                     nCellCenter = nCellCenter + nCellIt->second.bodies[i]->pos;
            //                 }
            //                 nCellCenter = nCellCenter * (1.0f / nCellIt->second.count);

            //                 Vector2 nScreenPos = worldToScreen(nCellCenter, scale, center);

            //                 Color lineColor = {255, 190, 152, (unsigned char)adaptiveAlpha};
            //                 rlColor4ub(lineColor.r, lineColor.g, lineColor.b, lineColor.a);
            //                 rlVertex2f(screenPos.x, screenPos.y);
            //                 rlVertex2f(nScreenPos.x, nScreenPos.y);

            //                 totalConnections++;
            //             }
            //         }
            //     }
            // }
            // else
            // {

            for (int i = 0; i < cell.count; i++)
            {
                const Body *body1 = cell.bodies[i];
                int connectionsCount = 0;

                // Removed size culling for source body
                Vector2 pos1 = worldToScreen(body1->pos, scale, center);

                // Relaxed screen bounds check
                if (!isInScreenBounds(body1->pos, scale, center, MAX_DISTANCE))
                    continue;

                for (int dx = -1; dx <= 1 && connectionsCount < adaptiveMaxConnections; dx++)
                {
                    for (int dy = -1; dy <= 1 && connectionsCount < adaptiveMaxConnections; dy++)
                    {
                        uint64_t nKey = getGridKey(Vec2{
                            body1->pos.x + dx * adaptiveDistance,
                            body1->pos.y + dy * adaptiveDistance});

                        const auto &nCellIt = spatialGrid.find(nKey);
                        if (nCellIt != spatialGrid.end())
                        {
                            const auto &nCell = nCellIt->second;
                            for (int j = 0; j < nCell.count && connectionsCount < adaptiveMaxConnections; j++)
                            {
                                const Body *body2 = nCell.bodies[j];
                                if (body1 != body2)
                                {
                                    // Removed size culling for target body
                                    float distSq = (body1->pos - body2->pos).mag_sq();
                                    if (distSq < adaptiveDistanceSq)
                                    {
                                        Vector2 pos2 = worldToScreen(body2->pos, scale, center);

                                        // Only check if target is in view
                                        if (!isInScreenBounds(body2->pos, scale, center, MAX_DISTANCE))
                                            continue;

                                        float alpha = std::max(0.0f, 1.0f - (std::sqrt(distSq) / adaptiveDistance));
                                        // Color lineColor = {255, 255, 255, (unsigned char)(alpha * adaptiveAlpha)};
                                        Color lineColor = {255, 0, 0, (unsigned char)(alpha * adaptiveAlpha)};

                                        rlColor4ub(lineColor.r, lineColor.g, lineColor.b, lineColor.a);
                                        rlVertex2f(pos1.x, pos1.y);
                                        rlVertex2f(pos2.x, pos2.y);

                                        connectionsCount++;
                                        totalConnections++;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    rlEnd();

    renderStats.lastFrameTime = static_cast<float>(GetTime()) - frameStart;
    renderStats.totalConnections = totalConnections;
}
struct QuadtreeRenderStats
{
    float lastFrameTime = 0.0f;
    float smoothedPerformanceRatio = 1.0f;
} quadtreeStats;


void drawQuadtreeNode(const Node &root, const std::vector<Node> &nodes, float scale, Vector2 center)
{
    if (!SHOW_QUADTREE)
        return;

    float frameStart = GetTime();

    // Pre-allocate vectors
    static std::vector<Rectangle> visibleQuads;
    static std::vector<Vector2> visiblePoints;
    static std::vector<const Node *> nodeStack;

    visibleQuads.clear();
    visiblePoints.clear();
    nodeStack.clear();

    visibleQuads.reserve(nodes.size());
    visiblePoints.reserve(nodes.size());
    nodeStack.reserve(nodes.size());

    nodeStack.push_back(&root);

    // Remove performance-based returns and adaptive calculations
    // Just collect visible nodes
    while (!nodeStack.empty())
    {
        const Node *current = nodeStack.back();
        nodeStack.pop_back();

        if (isInScreenBounds(current->data.quad.center, scale, center))
        {
            Vector2 screenPos = worldToScreen(current->data.quad.center, scale, center);
            float size = current->data.quad.size * scale;

            if (size >= 1.0f)
            { // Only add if visible
                visibleQuads.push_back({screenPos.x - size / 2,
                                        screenPos.y - size / 2,
                                        size,
                                        size});

                if (!current->is_branch())
                {
                    visiblePoints.push_back(screenPos);
                }
            }
        }

        if (current->is_branch())
        {
            for (int i = 0; i < 4; i++)
            {
                if (current->children + i < nodes.size())
                {
                    nodeStack.push_back(&nodes[current->children + i]);
                }
            }
        }
    }

    // Batch render quads
    rlBegin(RL_LINES);
    for (const auto &quad : visibleQuads)
    {
        Color quadColor = PERFORMANCE_MODE ? (Color){100, 100, 100, 50} : // Dimmer in performance mode
                              (Color){100, 100, 100, 100};                // Normal otherwise

        rlColor4ub(quadColor.r, quadColor.g, quadColor.b, quadColor.a);

        rlVertex2f(quad.x, quad.y);
        rlVertex2f(quad.x + quad.width, quad.y);
        rlVertex2f(quad.x + quad.width, quad.y);
        rlVertex2f(quad.x + quad.width, quad.y + quad.height);
        rlVertex2f(quad.x + quad.width, quad.y + quad.height);
        rlVertex2f(quad.x, quad.y + quad.height);
        rlVertex2f(quad.x, quad.y + quad.height);
        rlVertex2f(quad.x, quad.y);
    }
    rlEnd();

    quadtreeStats.lastFrameTime = static_cast<float>(GetTime()) - frameStart;
}

void drawBlackHole(const Body &body, float scale, Vector2 center)
{
    Vector2 screenPos = worldToScreen(body.pos, scale, center);
    float screenRadius = std::max(2.0f, body.radius * scale);

    for (int i = 4; i >= 0; i--)
    {
        float alpha = (1.0f - i / 4.0f) * 1.1f;
        DrawCircleGradient(
            screenPos.x, screenPos.y,
            screenRadius * (1.0f + i * 1.4f),
            (Color){255, 255, 237, (unsigned char)(alpha * 255)},
            (Color){0, 0, 0, 0});
    }

    float diskStartRadius = screenRadius * 2.1f;
    float diskEndRadius = screenRadius * 10.51f;
    // int segments = 10096;
    int segments = 5048;

    for (int i = 0; i < segments; i++)
    {
        float angle = (i * 300.0f / segments) * DEG2RAD;
        float nextAngle = ((i + 1) * 390.0f / segments) * DEG2RAD;

        float brightness = 1.4f + 1.0f * (10.5f + cosf(angle));
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
        {0, 0, 255, 255},     // deep blue (Hyper-giant Blue Stars)
        {100, 100, 255, 255}, // blue (Blue Stars)
        {173, 216, 230, 255}, // light blue (Blue-White Stars)
        {219, 233, 244, 255}, // bluish white (White Stars)
        {255, 255, 200, 255}, // light yellow (Transition to White)
        {255, 240, 150, 255}, // yellow (Sun-like stars)
        {255, 150, 50, 255},  // light orange (Transition to Yellow)
        {255, 100, 0, 255},   // deep orange (Orange Dwarfs)
        {255, 50, 0, 255},    // orange red (Red Dwarfs)
        {200, 0, 0, 255},     // deep red (Brown Dwarfs)

        // {200, 0, 0, 255},     // deep red (Brown Dwarfs)
        // {255, 50, 0, 255},    // orange red (Red Dwarfs)
        // {255, 100, 0, 255},   // deep orange (Orange Dwarfs)
        // {255, 150, 50, 255},  // light orange (Transition to Yellow)
        // {255, 240, 150, 255}, // yellow (Sun-like stars)
        // {255, 255, 200, 255}, // light yellow (Transition to White)
        // {219, 233, 244, 255}, // bluish white (White Stars)
        // {173, 216, 230, 255}, // light blue (Blue-White Stars)
        // {100, 100, 255, 255}, // blue (Blue Stars)
        // {0, 0, 255, 255},     // deep blue (Hyper-giant Blue Stars)
        {0, 0, 2, 1} // invisible/black (neutron stars or potentially extreme cases)
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
            SHARED_BODIES = simulation->bodies;
            SHARED_QUADTREE = simulation->quadtree.nodes;
        }

        std::this_thread::sleep_for(std::chrono::milliseconds(1));
        // std::unique_lock<std::mutex> lock(simulationMutex);
        // frameCondition.wait(lock, []
        //                     { return frameRendered; });
        // frameRendered = false;
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

    InitWindow(WIDTH, HEIGHT, "N-Body Simulation");
    GuiLoadStyleDefault();
    // SetTargetFPS(60);
    // SetTargetFPS(GetFPS() < 30 ? 60 : 0);  // 0 means uncapped

    RenderTexture2D circleTexture = LoadRenderTexture(32, 32);
    BeginTextureMode(circleTexture);
    ClearBackground(BLANK);
    DrawCircle(circleTexture.texture.width / 2, circleTexture.texture.height / 2, circleTexture.texture.width / 2, WHITE);
    EndTextureMode();

    auto simulation = std::make_shared<Simulation>();

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
            SHOW_QUADTREE = !SHOW_QUADTREE;
        }
        if (IsKeyPressed(KEY_C))
        {
            SHOW_CONNECTIONS = !SHOW_CONNECTIONS;
        }
        if (IsKeyPressed(KEY_V))
        {
            SHOW_BODIES = !SHOW_BODIES; // vanish bodies
        }
        Rectangle sliderBounds = {10, 160, 140, 20};
        float currentDt = SIMULATION_DT.load();
        float sliderDt = std::min(currentDt, SLIDER_MAX_DT); // clamp display value

        // Handle T/Y keys with clamping to 0.1 before going red
        if (IsKeyPressed(KEY_T))
        {
            float newDt = currentDt * 1.5f;
            if (currentDt < SLIDER_MAX_DT)
            {
                newDt = std::min(newDt, SLIDER_MAX_DT);
            }
            SIMULATION_DT.store(newDt);
        }
        if (IsKeyPressed(KEY_Y))
        {
            SIMULATION_DT.store(currentDt * 0.666f);
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
        if (IsKeyPressed(KEY_P))
        {
            PERFORMANCE_MODE = !PERFORMANCE_MODE;
        }

        BeginDrawing();
        ClearBackground(BLACK);
        BeginMode2D(camera);

        // draw bodies
        {
            if (SHOW_CONNECTIONS)
            {
                drawConnections(SHARED_BODIES, scale, center); // massive bottleneck performance lost is perporttotal to how long the simulation has been running
            }

            if (SHOW_QUADTREE && !SHARED_QUADTREE.empty())
            {
                drawQuadtreeNode(SHARED_QUADTREE[0], SHARED_QUADTREE, scale, center);
            }
            std::lock_guard<std::mutex> lock(UPDATE_LOCK);
            if (SHOW_BODIES)
            {

                if (PERFORMANCE_MODE)
                {
                    const Body *centralMass = nullptr;
                    float maxRadius = 0;

                    std::vector<const Body *> visibleBodies;
                    visibleBodies.reserve(SHARED_BODIES.size());

                    for (const auto &body : SHARED_BODIES)
                    {
                        if (body.radius > maxRadius)
                        {
                            maxRadius = body.radius;
                            centralMass = &body;
                        }

                        if (isInScreenBounds(body.pos, scale, center))
                        {
                            visibleBodies.push_back(&body);
                        }
                    }

                    if (!visibleBodies.empty())
                    {
                        for (const Body *body : visibleBodies)
                        {
                            Vector2 screenPos = worldToScreen(body->pos, scale, center);
                            float screenRadius = std::max(1.0f, body->radius * scale);

                            if (screenRadius < 1.0f)
                            {
                                continue;
                            }

                            Color color = (body == centralMass) ? RED : WHITE;
                            // DrawCircleV(screenPos, screenRadius, color);
                            DrawTexturePro(
                                circleTexture.texture,
                                (Rectangle){0, 0, (float)circleTexture.texture.width, (float)circleTexture.texture.height},
                                (Rectangle){screenPos.x - screenRadius, screenPos.y - screenRadius, screenRadius * 2, screenRadius * 2},
                                (Vector2){screenRadius, screenRadius},
                                0.0f,
                                color);
                        }
                    }
                }
                else
                {

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

                    for (const auto &body : SHARED_BODIES)
                    {

                        if (!isInScreenBounds(body.pos, scale, center))
                        {
                            continue;
                        }

                        Vector2 screenPos = worldToScreen(body.pos, scale, center);
                        float screenRadius = std::max(1.0f, body.radius * scale);

                        if (screenRadius < 1.0f)
                        {
                            continue;
                        }

                        // Render the body
                        // DrawCircleV(screenPos, screenRadius, getStarColorWithBrightness(body, 3.0f)); ->massive bottleneck
                        DrawTexturePro(
                            circleTexture.texture,
                            (Rectangle){0, 0, (float)circleTexture.texture.width, (float)circleTexture.texture.height},
                            (Rectangle){screenPos.x - screenRadius, screenPos.y - screenRadius, screenRadius * 2, screenRadius * 2},
                            (Vector2){screenRadius, screenRadius},
                            0.0f,
                            getStarColorWithBrightness(body, 3.0f));
                    }

                    if (centralMass)
                    {
                        drawBlackHole(*centralMass, scale, center);
                    }
                }
            }
        }

        EndMode2D();

        DrawFPS(10, 100);

        // UI goodies
        DrawText("Bodies: ", 10, 10, 20, WHITE);
        DrawText(TextFormat(" %zu", SHARED_BODIES.size()), 80, 10, 20, GetFPS() >= 30 ? GREEN : (GetFPS() >= 15 ? ORANGE : RED));
        DrawText(TextFormat("Scale: %.2f", (double)scale), 10, 35, 20, WHITE);
        DrawText(PAUSED ? "PAUSED" : "RUNNING", 10, 60, 20, PAUSED ? RED : GREEN);

        Color originalBaseColor = GetColor(GuiGetStyle(SLIDER, BASE_COLOR_NORMAL));
        Color originalPressedColor = GetColor(GuiGetStyle(SLIDER, BASE_COLOR_PRESSED));
        Color originalBackgroundColor = GetColor(GuiGetStyle(SLIDER, BACKGROUND_COLOR));
        Color originalBorderColor = GetColor(GuiGetStyle(SLIDER, BORDER_COLOR_NORMAL));
        Color originalLineColor = GetColor(GuiGetStyle(SLIDER, LINE_COLOR));

        int originalSliderWidth = GuiGetStyle(SLIDER, SLIDER_WIDTH);
        int originalSliderPadding = GuiGetStyle(SLIDER, SLIDER_PADDING);

        GuiSetStyle(SLIDER, SLIDER_WIDTH, 10);
        GuiSetStyle(SLIDER, SLIDER_PADDING, 5);

        Color transparentBlack = {0, 0, 0, 128}; // RGBA, where A = 128 (50% transparency)

        if (currentDt > SLIDER_MAX_DT)
        {
            GuiSetStyle(SLIDER, BASE_COLOR_NORMAL, ColorToInt(RED));
            GuiSetStyle(SLIDER, BASE_COLOR_PRESSED, ColorToInt(RED));

            GuiSetStyle(SLIDER, BACKGROUND_COLOR, ColorToInt(transparentBlack));
            GuiSetStyle(SLIDER, BORDER_COLOR_NORMAL, ColorToInt(transparentBlack));

            GuiSetStyle(SLIDER, LINE_COLOR, ColorToInt(transparentBlack));
        }
        else
        {
            GuiSetStyle(SLIDER, BASE_COLOR_NORMAL, ColorToInt(transparentBlack));
            GuiSetStyle(SLIDER, BASE_COLOR_PRESSED, ColorToInt(WHITE));
            GuiSetStyle(SLIDER, BACKGROUND_COLOR, ColorToInt(transparentBlack));
            GuiSetStyle(SLIDER, BORDER_COLOR_NORMAL, ColorToInt(transparentBlack));
            GuiSetStyle(SLIDER, LINE_COLOR, ColorToInt(transparentBlack));
        }

        // Draw a semi-transparent background for the slider
        // DrawRectangleRec((Rectangle){5, 155, 150, 30}, transparentBlack);
        // DrawRectangleLinesEx((Rectangle){5, 155, 150, 30}, 1, WHITE);

        // Draw slider with clamped display value
        if (GuiSliderBar(sliderBounds, "dt", TextFormat("%.3f", (double)currentDt),
                         &sliderDt, SLIDER_MIN_DT, SLIDER_MAX_DT))
        {
            SIMULATION_DT.store(sliderDt);
        }

        GuiSetStyle(SLIDER, BASE_COLOR_NORMAL, ColorToInt(originalBaseColor));
        GuiSetStyle(SLIDER, BASE_COLOR_PRESSED, ColorToInt(originalPressedColor));
        GuiSetStyle(SLIDER, BACKGROUND_COLOR, ColorToInt(originalBackgroundColor));
        GuiSetStyle(SLIDER, BORDER_COLOR_NORMAL, ColorToInt(originalBorderColor));
        GuiSetStyle(SLIDER, LINE_COLOR, ColorToInt(originalLineColor));
        GuiSetStyle(SLIDER, SLIDER_WIDTH, originalSliderWidth);
        GuiSetStyle(SLIDER, SLIDER_PADDING, originalSliderPadding);

        DrawText("Controls:", 10, HEIGHT - 100, 20, WHITE);
        DrawText("WASD: Pan | R/F: Zoom | Space: Pause | Q: Toggle Quadtree", 10, HEIGHT - 75, 20, WHITE);
        DrawText("C: Show Edges | V: Vanish | P: Perf Mode | T/Y: Change Time", 10, HEIGHT - 55, 20, WHITE);

        DrawText(TextFormat("dt: %.3f", (double)SIMULATION_DT.load()), 10, 130, 20, WHITE);

        if (!PERFORMANCE_MODE)
        {
            {
                std::lock_guard<std::mutex> lock(UPDATE_LOCK);
                SystemMetrics metrics = calculateMetrics(SHARED_BODIES);

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
        }
        EndDrawing();
        // {
        //     std::lock_guard<std::mutex> lock(simulationMutex);
        //     frameRendered = true;
        //     frameCondition.notify_one();
        // }
    }

    UnloadRenderTexture(circleTexture);
    CloseWindow();
    return 0;
}