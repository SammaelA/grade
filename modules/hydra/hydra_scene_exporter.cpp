#include "hydra_scene_exporter.h"
#include "HydraAPI/hydra_api/HydraAPI.h"
using pugi::xml_node;
#include "core/scene.h"
#include <unistd.h>
#include <signal.h>
#include <iostream>
#include <vector>
#include <GLFW/glfw3.h>
#include <errno.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
namespace h_inter
{
    void InfoCallBack(const wchar_t* message, const wchar_t* callerPlace, HR_SEVERITY_LEVEL a_level)
{
  if (a_level >= HR_SEVERITY_WARNING)
  {
    if (a_level == HR_SEVERITY_WARNING)
      std::wcerr << L"WARNING: " << callerPlace << L": " << message; // << std::endl;
    else      
      std::wcerr << L"ERROR  : " << callerPlace << L": " << message; // << std::endl;
  }
}

void destroy()
{
  std::cout << "call destroy() --> hrSceneLibraryClose()" << std::endl;
  hrSceneLibraryClose();
}

#ifdef WIN32
BOOL WINAPI HandlerExit(_In_ DWORD fdwControl)
{
  exit(0);
  return TRUE;
}
#else
bool destroyedBySig = false;
void sig_handler(int signo)
{
  if(destroyedBySig)
    return;
  switch(signo)
  {
    case SIGINT : std::cerr << "\nmain_app, SIGINT";      break;
    case SIGABRT: std::cerr << "\nmain_app, SIGABRT";     break;
    case SIGILL : std::cerr << "\nmain_app, SIGINT";      break;
    case SIGTERM: std::cerr << "\nmain_app, SIGILL";      break;
    case SIGSEGV: std::cerr << "\nmain_app, SIGSEGV";     break;
    case SIGFPE : std::cerr << "\nmain_app, SIGFPE";      break;
    default     : std::cerr << "\nmain_app, SIG_UNKNOWN"; break;
    break;
  }
  std::cerr << " --> hrSceneLibraryClose()" << std::endl;
  hrSceneLibraryClose();
  destroyedBySig = true;
}
#endif
};

namespace hydra
{
  void prepare_hydra_export_settings_block(const Block &in_settings, Block &export_settings)
  {
    std::string file_name = in_settings.get_string("save_filename", "selection/result");
    int image_count = in_settings.get_int("image_count", 1);
    float distance = in_settings.get_double("distance", 150);
    float height = in_settings.get_double("height", 50);
    glm::ivec2 image_size = in_settings.get_ivec2("image_size", {512, 512});
    int rays_per_pixel = in_settings.get_int("rays_per_pixel", 64);
    bool render_terrain = in_settings.get_bool("render_terrain", true);
    
    Block cameras;
    for (int i = 0; i < image_count; i++)
    {
      Block camera_blk;
      camera_blk.add_vec3("camera_look_at", glm::vec3(0, height, 0));
      camera_blk.add_vec3("camera_pos", glm::vec3(distance * cos(2 * PI * i / ((float)image_count)), height, distance * sin(2 * PI * i / ((float)image_count))));
      camera_blk.add_vec3("camera_up", glm::vec3(0, 1, 0));
      cameras.add_block("camera", &camera_blk);
    }
    export_settings.add_block("cameras", &cameras);
    export_settings.add_int("image_width", image_size.x);
    export_settings.add_int("image_height", image_size.y);
    export_settings.add_int("rays_per_pixel", rays_per_pixel);
    export_settings.add_bool("need_terrain", render_terrain);
    export_settings.add_bool("white_terrain", true);
    export_settings.add_string("demo_copy_dir", "saves/" + file_name);
  }

  bool export_scene(std::string directory, Scene &scene, Block &export_settings)
  {
    pid_t pid = fork();
    if (pid == -1)
    {
      // something REALLY bad happened
      perror("fork");
      exit(1);
    }
    else if (pid > 0)
    {
      // wait for child process and return if it exited normally
      int rv = -2;
      wait(&rv);
      int status = WEXITSTATUS(rv);
      logerr("hydra finished with status %d", status);
      return status == 0;
    }
    else
    {
      // child process. Does hydra stuff and dies

      hrInfoCallback(&h_inter::InfoCallBack);

      hrErrorCallerPlace(L"main"); // for debug needs only

      atexit(&h_inter::destroy);
#if defined WIN32
      SetConsoleCtrlHandler(&HandlerExit, TRUE); // if some one kill console :)
      wchar_t NPath[512];
      GetCurrentDirectoryW(512, NPath);
#ifdef NEED_DIR_CHANGE
      SetCurrentDirectoryW(L"../../main");
#endif
      std::wcout << L"[main]: curr_dir = " << NPath << std::endl;
#else
      std::string workingDir = "dependencies/Hydra/HydraAPI/main";
      if (chdir(workingDir.c_str()) != 0)
        std::cout << "chdir failed: " << workingDir.c_str() << std::endl;

      char cwd[1024];
      if (getcwd(cwd, sizeof(cwd)) != nullptr)
        std::cout << "[main]: curr_dir = " << cwd << std::endl;
      else
        std::cout << "getcwd() error" << std::endl;

      {
        struct sigaction sigIntHandler;
        sigIntHandler.sa_handler = h_inter::sig_handler;
        sigemptyset(&sigIntHandler.sa_mask);
        sigIntHandler.sa_flags = SA_RESETHAND;
        sigaction(SIGINT, &sigIntHandler, NULL);
        sigaction(SIGSTOP, &sigIntHandler, NULL);
        sigaction(SIGABRT, &sigIntHandler, NULL);
        sigaction(SIGILL, &sigIntHandler, NULL);
        sigaction(SIGTERM, &sigIntHandler, NULL);
        sigaction(SIGSEGV, &sigIntHandler, NULL);
        sigaction(SIGFPE, &sigIntHandler, NULL);
      }
#endif

      std::cout << "sizeof(size_t) = " << sizeof(size_t) << std::endl;

      try
      {
        export_internal(directory, scene, export_settings);
      }
      catch (std::runtime_error &e)
      {
        std::cout << "std::runtime_error: " << e.what() << std::endl;
      }
      catch (std::exception &e)
      {
        std::cout << "unknown exception" << e.what() << std::endl;
      }

      hrErrorCallerPlace(L"main"); // for debug needs only

      hrSceneLibraryClose();

      glfwTerminate();

      exit(0);
    }
    debug("Scene successfully exported to hydra_scene/%s\n", directory.c_str());
    return true;
  }
}