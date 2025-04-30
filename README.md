Fluid simulation engine in colaboration between Cochta and Chocolive 24 as a Bachelor main Project for SAE Institute

Olivier TODO:
- Faire une classe Renderer -> peut être mieux
- Faire un sample graphique raster et un rt
- Faire un dossier PhysicsEngine avec les libs Common, Physics et PhysicsSamples

Constantin TODO:
- Attendre Olivier

# Build prerequisites.
## For Falcor:
- Windows 10 version 20H2 (October 2020 Update) or newer, OS build revision .789 or newer
- Visual Studio 2022
- Windows 10 SDK (10.0.19041.0) for Windows 10, version 2004
- A GPU which supports DirectX Raytracing, such as the NVIDIA Titan V or GeForce RTX
- NVIDIA driver 466.11 or newer

## For building
CMake version 3.15 or newer.

# Build instructions.
1. Clone Nvidia's Falcor framework repository.
2. Add a subdirectory to the CMakeLists of falcor Samples (found in Source/Samples/) with this CMake command: 
**add_subdirectory(3DFluidSimulationEngine)**
3. Clone this repository inside the Source/Samples/ folder.
4. If you are working with Visual Studio 2022, you can setup a native Visual Studio solution by running setup_vs2022.bat after cloning this repository. The solution files are written to build/windows-vs2022 and the binary output is located in build/windows-vs2022/bin.
5. You can know run the visual studio solution and set the the 3DFluidSimulationEngine as the startup project.