#include "finiteDifferenceSolve.h"
#include "analyticalSolutions.h"
#include "SolvedElectrostaticSystem.h"
#include "UnsolvedElectrostaticSystem.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>  // Not needed?
#include <utility>
#include <sstream>
#include <unordered_map>
#include <cstdlib>
#include "sys/types.h"
#include "sys/sysinfo.h"
#include "sys/times.h"
#include "sys/vtimes.h"
#include "stdlib.h"
#include "stdio.h"

//parse line function 
int parseLine(char* line){
  int i =strlen(line);
  while (*line < '0' || *line > '9')line++;
  line[i-3] = '\0';
  i = atoi(line);
  return i;
}


//Returns RAM used by process in KB
int getRAMValue(){
  FILE* file = fopen("/proc/self/status","r");
  int result = -1;
  char line[128];
  while (fgets(line, 128, file) != NULL){
    if (strncmp(line, "VmRSS:", 6) == 0){
      result=parseLine(line);
      break;
    }
  }
  fclose(file);
  return result;
}

//Return CPU used by process
static clock_t lastCPU, lastSysCPU, lastUserCPU;
static int numProcessors;
    
    void init(){
        FILE* file;
        struct tms timeSample;
        char line[128];
    
        lastCPU = times(&timeSample);
        lastSysCPU = timeSample.tms_stime;
        lastUserCPU = timeSample.tms_utime;
    
        file = fopen("/proc/cpuinfo", "r");
        numProcessors = 0;
        while(fgets(line, 128, file) != NULL){
            if (strncmp(line, "processor", 9) == 0) numProcessors++;
        }
        fclose(file);
    }
    
    double getCurrentCPUValue(){
        struct tms timeSample;
        clock_t now;
        double percentCPU;
    
        now = times(&timeSample);
        if (now <= lastCPU || timeSample.tms_stime < lastSysCPU ||
            timeSample.tms_utime < lastUserCPU){
            percentCPU = -1.0;
        }
        else{
            percentCPU = (timeSample.tms_stime - lastSysCPU) +
                (timeSample.tms_utime - lastUserCPU);
            percentCPU /= (now - lastCPU);
            percentCPU /= numProcessors;
            percentCPU *= 100;
        }
        lastCPU = now;
        lastSysCPU = timeSample.tms_stime;
        lastUserCPU = timeSample.tms_utime;
    
        return percentCPU;
    }


void processLine(std::string &line, std::vector<std::string> &splitLine);

// A pair of electrostatic systems - one solved, one unsolved
struct systemPair {
    electrostatics::UnsolvedElectrostaticSystem unsolved;
    electrostatics::SolvedElectrostaticSystem solved;
    systemPair(int iMin, int iMax, int jMin, int jMax) :
        unsolved(iMin, iMax, jMin, jMax),
        solved(iMin, iMax, jMin, jMax) {
    }
};

/* TODO - WRITE THIS COMMENT
 */
int main(int argc, char* argv[]) {
    std::ifstream configFile;
    configFile.open(argv[1]);
    if(!configFile) {
        std::cerr << "Error Opening File\n";
        return(1);
    }

    // System pairs are indexed by names which the user enters
    std::unordered_map<std::string, systemPair> systems;
    std::string lastMade;   // Name of the last system pair made

    // Variables for reading a line
    std::string line;
    std::vector<std::string> splitLine;

    // Read in and process lines
    while(!configFile.eof()) {
        // Read a line and process it (split at spaces, convert to lowercase etc)
        std::getline(configFile, line);
        if(line == "") continue;             // Blank line
        splitLine.clear();
        processLine(line, splitLine);
        if(splitLine[0] == "#") continue;    // Comment line

        // For creating a new system
        else if(splitLine[0] == "new") {
            std::string name = splitLine[1];
            int iMin = std::stoi(splitLine[2]);
            int iMax = std::stoi(splitLine[3]);
            int jMin = std::stoi(splitLine[4]);
            int jMax = std::stoi(splitLine[5]);
            systems.emplace(name, systemPair(iMin, iMax, jMin, jMax));
            lastMade = name;
        }
            
        // For creating problem 1 analytical solution
        else if(splitLine[0] == "analytical1") {
            int iMin = systems.at(lastMade).solved.getIMin();
            int iMax = systems.at(lastMade).solved.getIMax();
            int jMin = systems.at(lastMade).solved.getJMin();
            int jMax = systems.at(lastMade).solved.getJMax();
            for(int i=iMin; i<= iMax; i++) {
                for(int j=jMin; j<=jMax; j++) {
                    double potential = electrostatics::analyticalProblem1(i, j, std::stod(splitLine[1]),
                            std::stod(splitLine[2]), std::stod(splitLine[3]), std::stod(splitLine[4]));
                    systems.at(lastMade).solved.setPotentialIJ(i, j, potential);
                }
            }
        }
        // For creating problem 2 analytical solution
        else if(splitLine[0] == "analytical2") {
            int iMin = systems.at(lastMade).solved.getIMin();
            int iMax = systems.at(lastMade).solved.getIMax();
            int jMin = systems.at(lastMade).solved.getJMin();
            int jMax = systems.at(lastMade).solved.getJMax();
            for(int i=iMin; i<= iMax; i++) {
                for(int j=jMin; j<=jMax; j++) {
                    double uniformField = electrostatics::uniformField(iMin, iMax, std::stod(splitLine[1]),
                            std::stod(splitLine[2]));
                    double potential = electrostatics::analyticalProblem2(i, j, std::stod(splitLine[3]), 
                            uniformField);
                    systems.at(lastMade).solved.setPotentialIJ(i, j, potential);
                }
            }
        }

        // For adding boundary conditions
        else if(splitLine[0] == "point") {
            systems.at(lastMade).unsolved.setBoundaryPoint(std::stoi(splitLine[1]), std::stoi(splitLine[2]),
                    std::stod(splitLine[3]));
        }
        else if(splitLine[0] == "ring") {
            systems.at(lastMade).unsolved.setBoundaryRing(std::stoi(splitLine[1]), std::stoi(splitLine[2]),
                    std::stod(splitLine[3]), std::stod(splitLine[4]));
        }
        else if(splitLine[0] == "circle") {
            systems.at(lastMade).unsolved.setBoundaryCircle(std::stoi(splitLine[1]), std::stoi(splitLine[2]),
                    std::stod(splitLine[3]), std::stod(splitLine[4]));
        }
        else if(splitLine[0] == "line") {
            systems.at(lastMade).unsolved.setBoundaryLine(std::stoi(splitLine[1]), std::stoi(splitLine[2]),
                    std::stoi(splitLine[3]), std::stoi(splitLine[4]), std::stod(splitLine[5]));
        }
        else if(splitLine[0] == "left") {
            systems.at(lastMade).unsolved.setLeftBoundary(std::stod(splitLine[1]));
        }
        else if(splitLine[0] == "right") {
            systems.at(lastMade).unsolved.setRightBoundary(std::stod(splitLine[1]));
        }
        else if(splitLine[0] == "top") {
            systems.at(lastMade).unsolved.setTopBoundary(std::stod(splitLine[1]));
        }
        else if(splitLine[0] == "bottom") {
            systems.at(lastMade).unsolved.setBottomBoundary(std::stod(splitLine[1]));
        }
        else if(splitLine[0] == "rectangle") {
            systems.at(lastMade).unsolved.setBoundaryRectangle(std::stoi(splitLine[1]), std::stoi(splitLine[2]),
                    std::stoi(splitLine[3]), std::stoi(splitLine[4]), std::stod(splitLine[5]));
        }

        // For solving - need to add options for different methods once they are properly implemented
        else if(splitLine[0] == "solve") {
            electrostatics::finiteDifferenceSolve(systems.at(lastMade).unsolved, systems.at(lastMade).solved);
        }

        // For comparisons and outputing results
        else if(splitLine[0] == "savesolution") {
            systems.at(lastMade).solved.saveFile(splitLine[1]);
        }
        // Would be better to implement a method for this with saving directly to a file in the future
        else if(splitLine[0] == "savecomparison") {
            int iMin = systems.at(splitLine[1]).solved.getIMin();
            int iMax = systems.at(splitLine[1]).solved.getIMax();
            int jMin = systems.at(splitLine[1]).solved.getJMin();
            int jMax = systems.at(splitLine[1]).solved.getJMax();
            electrostatics::SolvedElectrostaticSystem comparisonResult(iMin, iMax, jMin, jMax);
            systems.at(splitLine[1]).solved.compareTo(systems.at(splitLine[2]).solved, comparisonResult);
            comparisonResult.saveFile(splitLine[3]);
        }
        else if(splitLine[0] == "savefield") {
            systems.at(splitLine[1]).solved.saveFieldGNUPlot(splitLine[2]);
        }
    }
    configFile.close();
}


/* Process a line by splitting it at spaces, converting to lowercase, and saving 
 * each part in a vector.
 */
void processLine(std::string &line, std::vector<std::string> &splitLine) {
    std::transform(line.begin(), line.end(), line.begin(), ::tolower);
    std::istringstream lineStream(line);
    std::string word;
    while(std::getline(lineStream, word, ' ')) splitLine.push_back(word);
}
