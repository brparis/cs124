#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <stdint.h>
#include <vector>
#include <cmath>
#include <random>
#include <iostream>
#include <fstream>
#include <algorithm> // only for personal calculations (for table)
 
using namespace std;


/////////////////////////////////////////////////////////////////////////
// BEFORE SUBMITTING:                                                  //
// INPUT    ./randmst flag_of_choice numpoints numtrials dimension     //
// OUTPUT   average numpoints numtrials dimension                      //
/////////////////////////////////////////////////////////////////////////



//////////////////////
// helper functions //
//////////////////////
 
// Find the index of minimum element of float vector given a set of indices from which to check
// Remove the visited node (now part of MST)
// O(|V|)
int min_ind(vector<int> &indices, vector<float> &dists) {
    int ind = 0;                        // initialize to the first element in the indices vector in case none else exist and don't enter forloop
    float min = dists[indices[ind]];    // initialize minimum in case none else exist and don't enter for loop
    for(int i = 1; i < indices.size(); i++) {
        float curr = dists[indices[i]]; // only check answers of interest
        if(curr < min) {
            min = curr;
            ind = i;
        }
    }
    int retval = indices[ind];          // temporary storage of parent node (to be removed from "unvisited" list in next line)
    
    // remove now visited node (via faster swap-pop_back trick)
    iter_swap(indices.begin() + ind, indices.end() - 1);
    indices.pop_back();
 
    return retval;
}
 
 
// generate vector of random reals between 0 and 1 of specified input length
vector<float> rand01(int randlen) {

    /////////////////////////////
    // prep our rand generator //
    /////////////////////////////

    static random_device seeder;           // uses device to create random seed
    static mt19937 generator(seeder());    // Standard mersenne_twister_engine seeded with random_device seeder
    static uniform_real_distribution<double> dis(0.0, 1.0);    // distribution on which we randomly select

    vector<float> temprand;
    for(int i = 0; i < randlen; i++) {
        temprand.push_back(dis(generator));
    }
    return temprand;
}

// return a vector of euclidian distances (2D)
vector<float> euclido(vector<float> &xcoords, vector<float> &ycoords, vector<int> &unvisited, int currparent) {
    vector<float> distout;
    float x0 = xcoords[currparent];
    float y0 = ycoords[currparent];
    for(int i = 0; i < unvisited.size(); i++) {
        float ans = sqrt(pow((xcoords[unvisited[i]]-x0),2) + pow(ycoords[unvisited[i]]-y0, 2));
        distout.push_back(ans);
    }
    return distout;
}

// return a vector of euclidian distances (3D)
vector<float> euclido3d(vector<float> &xcoords, vector<float> &ycoords, vector<float> &zcoords, vector<int> &unvisited, int currparent) {
   vector<float> distout;
   float x0 = xcoords[currparent];
   float y0 = ycoords[currparent];
   float z0 = zcoords[currparent];
   for(int i = 0; i < unvisited.size(); i++) {
       float ans = sqrt(pow((xcoords[unvisited[i]]-x0),2) + pow(ycoords[unvisited[i]]-y0, 2)+pow((zcoords[unvisited[i]]-z0),2));
       distout.push_back(ans);
   }
   return distout;
}

// return a vector of euclidian distances (4D)
vector<float> euclido4d(vector<float> &xcoords, vector<float> &ycoords, vector<float> &zcoords, vector<float> &wcoords, vector<int> &unvisited, int currparent) {
   vector<float> distout;
   float x0 = xcoords[currparent];
   float y0 = ycoords[currparent];
   float z0 = zcoords[currparent];
   float w0 = wcoords[currparent];
   for(int i = 0; i < unvisited.size(); i++) {
       float ans = sqrt(pow((xcoords[unvisited[i]]-x0),2) + pow(ycoords[unvisited[i]]-y0, 2) + pow((zcoords[unvisited[i]]-z0),2) + pow((wcoords[unvisited[i]]-w0),2));
       distout.push_back(ans);
   }
   return distout;
}

 
int main(int argc,char* argv[])
{

    //////////////////////////
    // command line parsing //
    //////////////////////////
    
    int numpoints = 0;
    int numtrials = 0;
    int dimension = 0;
 
    // remember, Dhilan, argv[0] is "./randmst" //
    if(argc != 5) {
        printf("Too few or too many command line arguments passed. You passed %d arguments!\n\n", argc);
        return 1;
    }
    
    else if (argc == 5) {
        // convert arguments to integers //
        numpoints = atoi(argv[2]);
        numtrials = atoi(argv[3]);
        dimension = atoi(argv[4]);
    }

    // FOR TRIAL RUNS PRINT TO FILE
    // open a file in write/append mode
    ofstream outfile;
    outfile.open("trials.txt", ios::app);

    // STATS TO PRINT
    int trial_count = 0;
    float total_distance = 0;
    float max_on_trials = 0;
    double cpu_time;

    while (trial_count < numtrials) {
        trial_count++;
        clock_t start_time, end_time;


        start_time = clock();

        //////////////////////////
        // initialization stage //
        //////////////////////////
        
        vector<float> distance;         // tracks shortest distances
        vector<int> parent;             // parent from which to travel and achieve minimum distance
        vector<int> unvisited_nodes;    // not in MST yet
        int currparent = -1;            // used to update parent vector
        
        // let's create a vector of our unvisited nodes, set all of our distances to infinity, and parents to null
        for (int i = 0; i < numpoints; i++) {
            unvisited_nodes.push_back(i);   // useful indices for later
            distance.push_back(2);          // equivalent to infinity
            parent.push_back(-1);           // equivalent to NIL
        }

        // starting node
        distance[0] = 0;         // by definition, arbitrary starting node; dist := 0 and parent := NIL


        //////////////////
        //   printing   //
        //////////////////
        /*
        cout << "\n////////////////////////////////////\n        NOW WE'RE INITIALIZING        \n////////////////////////////////////\n";
        cout << "printing unvisited nodes vector of size: " << numpoints << "\n";
        for(int i = 0; i < unvisited_nodes.size(); i++) {
            cout << (float) unvisited_nodes[i] << " ";
        }
        cout << "\n \n";

        cout << "printing distance nodes vector of size: " << numpoints << "\n";
        for(int i = 0; i < distance.size(); i++) {
            cout << (float) distance[i] << " ";
        }
        cout << "\n \n";

        cout << "printing parent nodes vector of size: " << numpoints << "\n";
        for(int i = 0; i < parent.size(); i++) {
            cout << (float) parent[i] << " ";
        }
        cout << "\n////////////////////////////////////\n    NOW WE'RE DONE INITIALIZING    \n////////////////////////////////////\n";
        cout << "\n \n";  
        */

        if (dimension == 1) {
            /////////////////////////
            // adapted Prim's Algo //
            /////////////////////////

            while(unvisited_nodes.size() > 0) {
                
                // cout << "printing distance nodes vector of size: " << numpoints << "\n";
                // for(int i = 0; i < distance.size(); i++) {
                //     cout << (float) distance[i] << " ";
                // }
                // cout << "\n \n";



                currparent = min_ind(unvisited_nodes, distance);                      // find min distance (which we add to MST)

                // printf("current minimum unvisited node index is %i\n", currparent);

                // cout << "printing unvisited nodes vector of size: " << unvisited_nodes.size() << "\n";
                // for(int i = 0; i < unvisited_nodes.size(); i++) {
                //     cout << (float) unvisited_nodes[i] << " ";
                // }

                // cout << "\n////////////////////////////////////\n  One iteration, AH AH AH. \n ////////////////////////////////////\n";
                // cout << "\n \n";

                

                // generate edges from the node we're at now to all other nodes (aside from the ones already in MST)
                vector<float> newedges = rand01(unvisited_nodes.size());
                for(int i = 0; i < unvisited_nodes.size(); i++) {
                    int currchild = unvisited_nodes[i];              // provides relevant index
                    if(distance[currchild] > newedges[i]) {          // update if distance smaller  
                        distance[currchild] = newedges[i];
                        parent[currchild] = currparent;              // update parent, too
                    }
                }
            }
        }

        if (dimension == 2) {

            // initializing coordinates

            vector<float> xcoords = rand01(numpoints);
            vector<float> ycoords = rand01(numpoints);

            // cout << "printing (x, y) coordinates vector of size: " << numpoints << "\n";
            // for(int i = 0; i < xcoords.size(); i++) {
            //     cout << "(" << (float) xcoords[i] << ", " << (float) ycoords[i] << "); ";
            // }
            // cout << "\n \n";


            /////////////////////////
            // adapted Prim's Algo //
            /////////////////////////

            while(unvisited_nodes.size() > 0) {
                
                // cout << "printing distance nodes vector of size: " << numpoints << "\n";
                // for(int i = 0; i < distance.size(); i++) {
                //     cout << (float) distance[i] << " ";
                // }
                // cout << "\n \n";



                currparent = min_ind(unvisited_nodes, distance);                      // find min distance (which we add to MST)
                // printf("current minimum unvisited node index is %i\n", currparent);

                // search_n_destroy(unvisited_nodes, currparent);                        // "deletemin"




                // cout << "printing unvisited nodes vector of size: " << unvisited_nodes.size() << "\n";
                // for(int i = 0; i < unvisited_nodes.size(); i++) {
                //     cout << (float) unvisited_nodes[i] << " ";
                // }

                // cout << "\n////////////////////////////////////\n  One iteration, AH AH AH. \n ////////////////////////////////////\n";
                // cout << "\n \n";


                vector<float> newedges = euclido(xcoords, ycoords, unvisited_nodes, currparent);
                for(int i = 0; i < unvisited_nodes.size(); i++) {
                    int currchild = unvisited_nodes[i];              // provides relevant index
                    if(distance[currchild] > newedges[i]) {          // update if distance smaller  
                        distance[currchild] = newedges[i];
                        parent[currchild] = currparent;              // update parent, too
                    }
                }
            }      
        }

        if (dimension == 3) {

            // initializing coordinates

            vector<float> xcoords = rand01(numpoints);
            vector<float> ycoords = rand01(numpoints);
            vector<float> zcoords = rand01(numpoints);

            // cout << "printing (x, y, z) coordinates vector of size: " << numpoints << "\n";
            // for(int i = 0; i < xcoords.size(); i++) {
            //     cout << "(" << (float) xcoords[i] << ", " << (float) ycoords[i] << ", " << (float) zcoords[i] << "); ";
            // }
            // cout << "\n \n";


            /////////////////////////
            // adapted Prim's Algo //
            /////////////////////////

            while(unvisited_nodes.size() > 0) {
                
                // cout << "printing distance nodes vector of size: " << numpoints << "\n";
                // for(int i = 0; i < distance.size(); i++) {
                //     cout << (float) distance[i] << " ";
                // }
                // cout << "\n \n";



                currparent = min_ind(unvisited_nodes, distance);                      // find min distance (which we add to MST)
                // printf("current minimum unvisited node index is %i\n", currparent);

                // search_n_destroy(unvisited_nodes, currparent);                        // "deletemin"




                // cout << "printing unvisited nodes vector of size: " << unvisited_nodes.size() << "\n";
                // for(int i = 0; i < unvisited_nodes.size(); i++) {
                //     cout << (float) unvisited_nodes[i] << " ";
                // }

                // cout << "\n////////////////////////////////////\n  One iteration, AH AH AH. \n ////////////////////////////////////\n";
                // cout << "\n \n";


                vector<float> newedges = euclido3d(xcoords, ycoords, zcoords, unvisited_nodes, currparent);
                for(int i = 0; i < unvisited_nodes.size(); i++) {
                    int currchild = unvisited_nodes[i];              // provides relevant index
                    if(distance[currchild] > newedges[i]) {          // update if distance smaller  
                        distance[currchild] = newedges[i];
                        parent[currchild] = currparent;              // update parent, too
                    }
                }
            }      
        }

        if (dimension == 4) {

            // initializing coordinates

            vector<float> xcoords = rand01(numpoints);
            vector<float> ycoords = rand01(numpoints);
            vector<float> zcoords = rand01(numpoints);
            vector<float> wcoords = rand01(numpoints);

            // cout << "printing (x, y, z, w) coordinates vector of size: " << numpoints << "\n";
            // for(int i = 0; i < xcoords.size(); i++) {
            //     cout << "(" << (float) xcoords[i] << ", " << (float) ycoords[i] << ", " << (float) zcoords[i] << ", " << (float) wcoords[i] << "); ";
            // }
            // cout << "\n \n";


            /////////////////////////
            // adapted Prim's Algo //
            /////////////////////////

            while(unvisited_nodes.size() > 0) {
                
                // cout << "printing distance nodes vector of size: " << numpoints << "\n";
                // for(int i = 0; i < distance.size(); i++) {
                //     cout << (float) distance[i] << " ";
                // }
                // cout << "\n \n";



                currparent = min_ind(unvisited_nodes, distance);                      // find min distance (which we add to MST)
                // printf("current minimum unvisited node index is %i\n", currparent);

                // search_n_destroy(unvisited_nodes, currparent);                        // "deletemin"




                // cout << "printing unvisited nodes vector of size: " << unvisited_nodes.size() << "\n";
                // for(int i = 0; i < unvisited_nodes.size(); i++) {
                //     cout << (float) unvisited_nodes[i] << " ";
                // }

                // cout << "\n////////////////////////////////////\n  One iteration, AH AH AH. \n ////////////////////////////////////\n";
                // cout << "\n \n";


                vector<float> newedges = euclido4d(xcoords, ycoords, zcoords, wcoords, unvisited_nodes, currparent);
                for(int i = 0; i < unvisited_nodes.size(); i++) {
                    int currchild = unvisited_nodes[i];              // provides relevant index
                    if(distance[currchild] > newedges[i]) {          // update if distance smaller  
                        distance[currchild] = newedges[i];
                        parent[currchild] = currparent;              // update parent, too
                    }
                }
            }      
        }

        end_time = clock();
        cpu_time += ( (double) (end_time - start_time) ) / (double) CLOCKS_PER_SEC;

        // collecting total distance
        for (int i = 0; i < distance.size(); i++) {
            total_distance += distance[i];
        }
        
        // updating max element
        if (*max_element(distance.begin(), distance.end()) > max_on_trials) {
            max_on_trials = *max_element(distance.begin(), distance.end());
        }

    }


    outfile << "///////////////////////////////////////////////////////////////////////////////////////" << endl;
    outfile << "CPU Time (Avg), Max Dist, Total Distance (Avg), Numpoints, Numtrials, Dimension" << endl;
    outfile << (double) cpu_time / (double) numtrials << endl;
    outfile << max_on_trials <<  endl;
    outfile << total_distance / (float) numtrials << endl;
    outfile << numpoints << endl;
    outfile << numtrials << endl;
    outfile << dimension <<endl;
    outfile << "///////////////////////////////////////////////////////////////////////////////////////" << endl;
    outfile.close();

    cout << total_distance / (float) numtrials << " " << numpoints << " " << numtrials << " " << dimension << endl;
    return 0;
}
