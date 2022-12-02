#include "AMDiS.h"
#include "PhaseFieldProblemSimple_divisionless_tensionless.h" 		
#include "GlobalProblem.h"
#include "Coupling_divisionless_tensionless_phase.h"


#include "AdaptInstationarySeq.h"

using namespace AMDiS;

#ifndef DIM
#define DIM 2
#endif

int main(int argc, char** argv)
{
  AMDiS::init(argc, argv);
  mpi14::Environment env(argc, argv);
  mpi14::Communicator world;

  Timer t;
  std::srand(std::time(0) + world.rank());
  
  int num_cells = 5;
  Parameters::get("number of cells", num_cells);
  //TEST_EXIT(world.size() == num_cells)("ntasks != ncells\n");

  std::string setup = "";
  Parameters::get("main->setup", setup);

  // Firstly catch wrong inputs of cell numbers
  if (setup.compare("singular") == 0)
    TEST_EXIT(num_cells == 1)("Singular only works for 1 cell.");

  if (setup.compare("binary_collision") == 0)
    TEST_EXIT(num_cells == 2)("Binary collision only works for 2 cells.");
  if (setup.compare("binary_Lea") == 0)
    TEST_EXIT(num_cells == 2)("binary_Lea only works for 2 cells.");
  if (setup.compare("three_collision") == 0)
    TEST_EXIT(num_cells == 3)("three_collision only works for 3 cells.");
  if (setup.compare("t1") == 0){
    double n = std::floor(std::sqrt(num_cells));
    TEST_EXIT(std::fmod(n,2.0) == 0)("Only squares of even numbers allowed!");
    TEST_EXIT(std::abs(std::sqrt(num_cells) - n ) < 1e-5)("Not a square number!");
  }

  base_problems::PhaseFieldProblem phaseProb{world};
  base_problems::GlobalProblem globalProb{"interaction",world};

  base_problems::PhaseFieldGlobal phaseGlobal{phaseProb, globalProb};
  phaseGlobal.initialize(INIT_ALL);

  AdaptInfo adaptInfo("adapt", phaseGlobal.getNumComponents());
  AdaptInstationarySeq adapt("adapt", phaseGlobal, adaptInfo, phaseGlobal, adaptInfo);
  
  phaseGlobal.initTimeInterface();
  adapt.adapt();

  MSG("time (total) = %e\n", t.elapsed());

  AMDiS::finalize();
}
