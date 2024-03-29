#include "JoojakApp.h"
#include "MooseInit.h"
#include "Moose.h"
#include "MooseApp.h"
#include "AppFactory.h"

// Create a performance log
PerfLog Moose::perf_log("Joojak");

/*
     _             _       _
    | | ___   ___ (_) __ _| | __
 _  | |/ _ \ / _ \| |/ _` | |/ /
| |_| | (_) | (_) | | (_| |   <
 \___/ \___/ \___// |\__,_|_|\_\
                |__/
*/

int main(int argc, char *argv[])
{
  // Initialize MPI, solvers and MOOSE
  MooseInit init(argc, argv);

  // Register this application's MooseApp and any it depends on
  JoojakApp::registerApps();

  // This creates dynamic memory that we're responsible for deleting
  MooseApp * app = AppFactory::createApp("JoojakApp", argc, argv);

//  app->legacyUoInitializationDefault() = false;
//  app->legacyUoAuxComputationDefault() = false;

  // Execute the application
  app->run();

  // Free up the memory we created earlier
  delete app;

  return 0;
}
