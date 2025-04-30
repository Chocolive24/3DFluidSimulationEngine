#include "FluidApplication.h"

int runMain(int argc, char** argv)
{
    SampleAppConfig config;
    config.windowDesc.title = "3D Fluid Simulation Engine";
    config.windowDesc.resizableWindow = true;

    FluidApplication helloDXR(config);
    return helloDXR.run();
}

int main(int argc, char** argv)
{
    return catchAndReportAllExceptions([&]() { return runMain(argc, argv); });
}
