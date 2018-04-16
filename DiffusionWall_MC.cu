/*Raul P. Pelaez 2018. 
  An UAMMD[1] Monte Carlo simulation of Debye type blobs confined near a wall with a Debye repulsion.
  It uses the MonteCarloNVT_Anderson module by Pablo Ibañez[2].

  ### USAGE
  
  ## HOW TO COMPILE THIS FILE

    Compile with nvcc, you need to tell it where UAMMD and HydroGrid[3] are (HydroGrid has to be compiled already) and the flag --expt-relaxed-constexpr is necessary.
    You will need to link HydroGrid in some way. I usually pass -lCallHydroGrid and then copy libCallHydroGrid.so to the same folder as the executable, it ight be different for your system.
    Refer to [5] for more compilation options and details. You should have recieved a Makefile along this file though.

    
  ## HOW TO RUN A SIMULATION
  
    This code needs the following files to be in the same folder in order to run:
    
      -data.main            -> A file containing the simulation parameters, you should have received an example data.main along this file.
      -hydroGridOptions.nml -> A HydroGrid configuration file. See [3] for more info.
      
    After compiling, just execute the binary. You can pass the --device X flag to select a certain GPU.
    
  ### SYSTEM DETAILS 

    Particles interact via a repulsive Debye potential.
  
    Particles interact with two walls in the xy plane.  
    The walls are:
           1- At Lz/2 (infinite wall)
  	   2- At -Lz/2(infinite wall + Debye repulsion)
    
  
    The system is technically unbounded in the z direction. The two inifinite zwalls avoid periodic artifacts like particles seeing each other through the floor/ceiling of the box. The box is increased by blob_cutOff in the z direction to avoid particles near one of the walls interact with the other due to PBC.
  
  ### PARAMETERS AND INITIAL CONFIGURATION

    All parameters are read from a data.main file. You should have recieved an example reference one with this file.
    You can refer to the function readParameters to see a list of available parameters.
    
    The initial configuration of particles will be random if a coordinates file is not provided in data.main
    The format of the coordinates file provided in the data.main is:
    ----
    Nparticles
    x y z
    ...
    ----
    A fourth column might be included in the file that will be interpreted as species by HydroGrid. If only three columns are available all particle are assumed to have species 0.

  ### HYDROGRID ANALYSIS
  
    This code is able to call HydroGrid[3] to process the simulation results on the fly.
    
    HydroGrid is called every samplefreq to update itself according to hydroGridOptions.nml (which must be present at the current folder) and write to disk. You should have received an example .nml with this file, see [3] for more information about HydroGrid and its input.
  
    HydroGrid is called through the UAMMD wrapper class HydroGrid in utils/HydroGrid.cuh, you can see more details there.
  
  ### OUTPUT

    Information about the current simulation will be output at runtime through stdout and stderr,
        including all parameters used and some useful information. 
  
    Besides HydroGrid output, this code will output the particle positions every savefreq with the following format:
    ----
    #time1
    x y z radius color
    ...
    #time2
    ...
    ----
    The positions will not be folded to the primary simulation box.
  
    You can visualize the particle results with superpunto[4]. HG will output whatever was asked in the .nml.
  

   ### References: 
   [1] https://github.com/RaulPPelaez/UAMMD/wiki
   [2] https://github.com/PabloIbannez/UAMMD/blob/master/src/Integrator/MonteCarlo/MonteCarloNVT_Anderson.cuh
   [3] https://github.com/stochasticHydroTools/HydroGrid
   [4] https://github.com/RaulPPelaez/superpunto
   [5] https://github.com/RaulPPelaez/UAMMD/wiki/Compiling-UAMMD
*/

//This include contains the basic needs for an uammd project
#include"uammd.cuh"
//The rest can be included depending on the used modules
#include"Integrator/MonteCarlo/MonteCarloNVT_Anderson.cuh"
#include"Interactor/Potential/Potential.cuh"
#include"utils/HydroGrid.cuh"
#include"utils/InitialConditions.cuh"
#include<fstream>
#include<sstream>
/*************  PARAMETERS AND DEFAULT VALUES *******************/
//All these parameters will be override if found in the data.main

int numberofparticles;
double3 boxSize;

int numsteps;   //Total number of simulation steps
int savefreq;   //The code will output every savefreq steps
int numstepsRelaxation;  //The code will run these steps as thermalization steps

bool loadParticles = false; //true will load initial configuration from "coordinates" file
std::string coordinates;    //Name of the file with the initial configuration
//if colors==true, color will be read from the fourth column in the input file and HydroGrid will interpret it as species.
//If false the coordinates file has to contain only three columns
bool useColors = false;

//Debye potential parameters for the blob-blob interaction, see Debye for more info
//These are the default values that will be override if provided in the data.main
double blob_repulsion_strength = 0.0165677856; 
double blob_radius = 0.656;
double blob_debye_length = 0.0656;
double blob_cutOff=1.968;

//Intensity of the gravity
double gravity = 0.0869477388288;

//Debye potential parameters of the blob-wall interaction
//These are the default values that will be override if provided in the data.main
double wall_repulsion_strength = 0.482995388566091;
double wall_blob_radius = 0.656;
double wall_debye_length = 2.624;

//System temperature
double temperature;

//desired acceptance ratio for the MCMC algorithm
double acceptanceRatio=0.5;
//Starting value of the random displacement of the particles for the MCMC algorithm (it will be autotuned to achieve the desired acceptance ratio)
double initialDx;

//Number of cells for the HydroGrid analysis
int3 cellsHydroGrid = make_int3(128, 128, 1);
//Steps to update HG with the current simulation state
int samplefreqHydroGrid = -1;
//All the simulation output files will be prepended by this
std::string outputName = "run";

/**********************************************************************************************/

//This function will parse a data.main file, reading all the above variables.
void readParameters(std::string datamain);

using namespace uammd;

//Interparticle potential. Force is not necessary for MC
//See how it is used below inside main, look for DebyePot
struct Debye{
  struct InputPairParameters{
    real cutOff, blob_radius, debye_length, repulsion_strength;
  };      
  using PairParameters = InputPairParameters;
  
  //Return F/r
  static inline __host__ __device__ real force(const real &r2, const PairParameters &params){
    if(r2 >= params.cutOff*params.cutOff) return 0;
    const real r = sqrt(r2);    
    real fmod;
    if(r<=real(2.0)*params.blob_radius)
      fmod = -params.repulsion_strength/(params.debye_length*r);
    else
      fmod = -params.repulsion_strength*exp(-(r-real(2.0)*params.blob_radius)/params.debye_length)/(r*params.debye_length);
    
    return fmod;  
  }
      
  static inline __host__ __device__ real energy(const real &r2, const PairParameters &params){
    const real r = sqrt(r2);    
    real E;
    if(r<=real(2.0)*params.blob_radius)
      E = params.repulsion_strength*(real(1.0) + (real(2.0)*params.blob_radius -r)/params.debye_length);
    else
      E = params.repulsion_strength*exp(-(r-real(2.0)*params.blob_radius)/params.debye_length);

    return E;      
  }



  static inline __host__ PairParameters processPairParameters(InputPairParameters in_par){return in_par;}
    
};
//DebyePot is a kind of radial potential
using DebyePot = Potential::Radial<Debye>;

//External potential
struct Wall{
  real gravity;
  real zwall_top, zwall_bottom;
  Debye::PairParameters par;
  Wall(real gravity, real zwall_bottom, real zwall_top, Debye::PairParameters par):
    gravity(gravity), zwall_top(zwall_top),zwall_bottom(zwall_bottom), par(par){}
    
  //Force is not needed for MCMC
  // inline __device__ real3 force(const real4 &pos){
  //   real3 f= make_real3(0, 0, -k);

  //   real rwall = zwall_top - pos.z;
  //   f += Debye::force(rwall*rwall, par)*make_real3(0, 0, rwall);
  //   rwall = zwall_bottom - pos.z;
  //   f += Debye::force(rwall*rwall, par)*make_real3(0, 0, rwall);
  //   return f;
    
  // }

  inline __device__ real energy(const real4 &pos){
        
    if(pos.z<zwall_bottom || pos.z>zwall_top) return INFINITY;
    real rwall = pos.z-zwall_bottom;
    real E = gravity*rwall;

    E += Debye::energy(rwall*rwall, par);
    return E;
  }
};


using namespace std;

int main(int argc, char *argv[]){
//Process data.main
  std::string inputFile="data.main";
  readParameters(inputFile);
  
  real lz = boxSize.z;

  real zwall_top = lz*0.5; //Location of the bottom wall
  real zwall_bottom = -lz*0.5;  //Location of the top wall


  //UAMMD System entity holds information about the GPU and tools to interact with the computer itself (such as a loging system). All modules need a System to work on.  
  auto sys = make_shared<System>();

  //Modules will ask System when they need a random number (i.e for seeding the GPU RNG).
  ullint seed = 0xf31337Bada55D00dULL^time(NULL);
  sys->rng().setSeed(seed);

  //ParticleData stores all the needed particle properties the simulation will need (pos, vel...).
  auto pd = make_shared<ParticleData>(numberofparticles, sys);

  boxSize.z += blob_cutOff;
  Box box(make_real3(boxSize));
  //Initial positions
  {
    //Ask pd for a property like so:
    auto pos = pd->getPos(access::location::cpu, access::mode::write);
    
    //Load particles from a file
    if(loadParticles){
      sys->log<System::MESSAGE>("Reading initial positions from %s", coordinates.c_str());
      //Count columns
      int ncolumns = -1;
      {
	ifstream in(coordinates);
			
	string line;
	getline(in, line);
	getline(in, line);
	istringstream iss(line);
	ncolumns = (vector<string>{istream_iterator<string>{iss},
	      istream_iterator<string>{}}).size();

      }
      useColors = ncolumns==4;
      ifstream in(coordinates);
      fori(0, numberofparticles){
	in>>pos.raw()[i].x>>pos.raw()[i].y>>pos.raw()[i].z;
	//Type of particle is stored in .w
	pos.raw()[i].w = 0;
	if(useColors)
	  in>>pos.raw()[i].w;
      }
      if(useColors)
	sys->log<System::MESSAGE>("Using fourth column of %s as species", coordinates.c_str());
      else
	sys->log<System::MESSAGE>("No species column detected in %s", coordinates.c_str());
	    

    }
    //Start randomly between the two walls
    else{
      sys->log<System::MESSAGE>("Generating random initial positions");
      fori(0, numberofparticles){
	pos.raw()[i] = make_real4(boxSize*sys->rng().uniform3(-0.5, 0.5) , 0.0);
	pos.raw()[i].z = sys->rng().uniform(zwall_bottom, zwall_top);
	pos.raw()[i].w = sys->rng().next32()%2;
      }
    }
  }

  sys->log<System::MESSAGE>("%d particles created.", numberofparticles);

  //Most modules work with particle groups, if no criteria is provided the group contains all particles, like in this case.
  auto pg = make_shared<ParticleGroup>(pd, sys, "All");

  //Blob blob interaction
  auto pot = make_shared<DebyePot>(sys);
 
  Debye::InputPairParameters blobPar;
  blobPar.repulsion_strength = blob_repulsion_strength;
  blobPar.blob_radius = blob_radius;
  blobPar.debye_length = blob_debye_length;
  blobPar.cutOff= blob_cutOff;

  //Each Potential describes the pair interactions with certain parameters.
  //The first two parameters are the particle colors, if only 0,0 is passed the parameters will be used for all colors.
  //The needed ones are in InputPairParameters inside each potential, in this case:  
  pot->setPotParameters(0, 0, blobPar);
    
  //External potential (gravity + Wall)
  real g = gravity;
  
  Debye::InputPairParameters wallRep;
  wallRep.repulsion_strength = wall_repulsion_strength;
  wallRep.blob_radius = wall_blob_radius;
  wallRep.debye_length = wall_debye_length;
  wallRep.cutOff=boxSize.x*2.0;
  
  auto wall = make_shared<Wall>(g, zwall_bottom, zwall_top, wallRep);
  
  //Monte Carlo parameters
  using MC = MC_NVT<DebyePot,Wall>;
  MC::Parameters par;
  par.box = box;         //Simulation box
  par.kT = temperature;  //temperature
  par.attempsNumPerCell = 10;  //Number of movement tries in each cell
  par.shiftValue = blob_cutOff*1.2;   //Grid stride length
  par.posVar = initialDx;              //Starting particle jump
  par.thermSteps = numstepsRelaxation;       //Number of self optimization steps
  par.accRatio = acceptanceRatio;         //Desired acceptance ratio
  par.accFactor = 1.11;       //Change in jump every optimization step
  par.checkAccRate = 10;       //Check acceptance every X during thermalization

  //Create the MonteCarlo UAMMD integrator object.
  auto mc = make_shared<MC>(pd, pg, sys, pot, wall, par);
      

   
  //You can issue a logging event like this, a wide variety of log levels exists (see System.cuh).
  //A maximum log level is set in System.cuh, every logging event with a level superior to the max will result in
  // absolutely no overhead, so dont be afraid to write System::DEBUGX log calls.
  sys->log<System::MESSAGE>("RUNNING!!!");

  //This can measure time
  Timer tim;
  
  ofstream out(outputName + string(".particle.pos")); //File to output particle positions

  //HydroGrid configuration
  HydroGrid::Parameters hgpar;
  hgpar.box = box;               //Simulation box
  hgpar.cellDim = cellsHydroGrid;//cells to perform HG analysis
  hgpar.dt = 1.0;                //Time between update calls
  hgpar.outputName = outputName; //Name prefix of HG output
  hgpar.useColors = true;        //Interpret pos.w as species

  HydroGrid hg(pd, sys, hgpar);

  //Initialize HydroGrid
  hg.init();
  
  //Thermalization. the first par.thermSteps calls to mc->forwardTime() will be perform as thermalization steps
  forj(0,par.thermSteps){mc->forwardTime();}
  
  tim.tic();
  //Run the simulation
  forj(0,numsteps){
    //This will instruct the integrator to take the simulation to the next time step,
    //For the MCMC module this means doing a certain number of trials[2].
    mc->forwardTime();

    //update HG
    if(samplefreqHydroGrid>0 && j%samplefreqHydroGrid==0) {hg.update(j);} 
    
    //Write results
    if(j%savefreq==0 && savefreq > 0)
    {
      sys->log<System::MESSAGE>("[System] Writing to disk...");

      //HG write current results
      if(samplefreqHydroGrid > 0){hg.write(j);}
      
      auto pos = pd->getPos(access::location::cpu, access::mode::read);
      //ParticleData has the right of changing the particles order at any time.
      //This allows to always recover the original order of the particles.
      const int * sortedIndex = pd->getIdOrderedIndices(access::location::cpu);
      
      out<<"#\n";
      real3 p;
      fori(0, numberofparticles){
	real4 pc = pos.raw()[sortedIndex[i]];
	//p = box.apply_pbc(make_real3(pc));
	p = make_real3(pc);
	int type = pc.w;
	out<<p<<" "<<0.5<<" "<<type<<"\n";
      }
      out<<flush;
    }
  }
  
  auto totalTime = tim.toc();
  sys->log<System::MESSAGE>("mean FPS: %.2f", numsteps/totalTime);
  //sys->finish() will ensure a smooth termination of any UAMMD module.
  sys->finish();

  return 0;
}



//Read parameters from a file, you can include any new parameter here
void readParameters(std::string datamain){
  std::ifstream in(datamain);

  std::string word, equal;
  while(in>>word){
    if(word.compare("numberofparticles")==0)        in>>equal>>numberofparticles;
    else if(word.compare("cutoff")==0)              in>>equal>>blob_cutOff;
    else if(word.compare("temperature")==0)         in>>equal>>temperature;
    else if(word.compare("numstepsRelaxation")==0)  in>>equal>>numstepsRelaxation;
    else if(word.compare("cellsHydroGrid")==0)      in>>equal>>cellsHydroGrid.x>>cellsHydroGrid.y>>cellsHydroGrid.z;
    else if(word.compare("numsteps")==0)            in>>equal>>numsteps;
    else if(word.compare("initialDx")==0)           in>>equal>>initialDx;
    else if(word.compare("acceptanceRatio")==0)     in>>equal>>acceptanceRatio;
    else if(word.compare("samplefreq")==0)          in>>equal>>samplefreqHydroGrid;
    else if(word.compare("savefreq")==0)            in>>equal>>savefreq;
    else if(word.compare("loadParticles")==0)       in>>equal>>loadParticles;
    else if(word.compare("lx")==0)                  in>>equal>>boxSize.x;
    else if(word.compare("ly")==0)                  in>>equal>>boxSize.y;
    else if(word.compare("lz")==0)                  in>>equal>>boxSize.z;
    else if(word.compare("coordinates")==0)         in>>equal>>coordinates;
    else if(word.compare("blob_repulsion_strength")==0) in>>equal>>blob_repulsion_strength;
    else if(word.compare("blob_radius")==0)	       in>>equal>>blob_radius;
    else if(word.compare("blob_debye_length")==0)       in>>equal>>blob_debye_length;
    else if(word.compare("blob_cutOff")==0)	       in>>equal>>blob_cutOff;
    else if(word.compare("gravity")==0)		       in>>equal>>gravity;
    else if(word.compare("wall_repulsion_strength")==0) in>>equal>>wall_repulsion_strength;
    else if(word.compare("wall_blob_radius")==0)	       in>>equal>>wall_blob_radius;
    else if(word.compare("wall_debye_length")==0)       in>>equal>>wall_debye_length;
    else if(word.compare("outputname")==0)          in>>equal>>outputName;
    else if(word.compare("initParticles")==0)       {getline(in,equal);}//This options does nothing
    else if(word.compare("wallX")==0)               {getline(in,equal);}//This options does nothing
    else if(word.compare("wallY")==0)               {getline(in,equal);}//This options does nothing
    else if(word.compare("wallZ")==0)               {getline(in,equal);}//This options does nothing
    else if(word.find_first_of("!#/&")!=std::string::npos) {getline(in,equal);} //Ignore comments
    else{
      cerr<<"ERROR: Unrecognized word in data.main!!"<<endl;
      cerr<<"Last word readed: "<<word<<endl;
      exit(1);
    }
  }
  in.close();
  
  //Copy input file to output 
  ofstream(outputName+"."+datamain)<<ifstream(datamain).rdbuf();
    
}
