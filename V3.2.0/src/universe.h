#ifndef CAC_UNIVERSE_H
#define CAC_UNIVERSE_H

#include <stdio.h>
#include <vector>
#include <string>
#include "pointers.h"

namespace CAC_NS {

// used for unordered_map with std::vector<int> as key
template <typename T>
struct VectorHasher {
  std::size_t operator()(const std::vector<T> &vec) const {
    std::size_t seed = vec.size();
    for (const T &i : vec) {
      // Combine hash values for each element
      seed ^= std::hash<T>()(i) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }
    return seed;
  }
};

// This struct is ALWAYS size 2, making it a fixed-size key.
struct Vector2D {
    int data[2];

    // Constructor to initialize the two elements
    Vector2D(int a, int b) {
        data[0] = a;
        data[1] = b;
    }

    // Required for map equality check
    bool operator==(const Vector2D& other) const {
        return data[0] == other.data[0] && data[1] == other.data[1];
    }
    
    // Default constructor
    Vector2D() : data{0, 0} {} 
};

// used for unordered_map with Vector2D as key
struct Vector2DHasher {
    size_t operator()(const Vector2D& key) const {
        size_t h1 = std::hash<int>{}(key.data[0]);
        size_t h2 = std::hash<int>{}(key.data[1]);
        
        // Hash combine logic
        return h1 ^ (h2 + 0x9e3779b9 + (h1 << 6) + (h1 >> 2));
    }
};

class Universe : protected Pointers {
 public:
  const char *version;    // CAC version string

  MPI_Comm uworld;        // communicator for entire universe
  int me, nprocs;          // my place in universe

  FILE *uscreen;          // universe screen output
  FILE *ulogfile;         // universe logfile

  int existflag;          // 1 if universe exists due to -partition flag
  int nworlds;            // # of worlds in universe
  int iworld;             // which world I am in
  int *procs_per_world;   // # of procs in each world
  int *root_proc;         // root proc in each world

  MPI_Comm uorig;         // original communicator passed to LAMMPS instance
  int *uni2orig;          // proc I in universe uworld is
                          // proc uni2orig[I] in original communicator

  Universe(class CAC *, MPI_Comm);
  ~Universe();
  void reorder(char *, char *);
  void add_world(char *);
  int consistent();

  // frequently used universal functions
  
  void bounds(const char *, int, char*, int, int &, int &, int nmin=1);
  void boundsbig(const char *, int, char *, bigint, bigint &, bigint &, bigint nmin=1);
  double numeric(const char *, int, char *);
  int inumeric(const char *, int, char *);
  bigint bnumeric(const char *, int, char *);
  tagint tnumeric(const char *, int, char *);
  int count_words(const char *);
  int next_prime(int n);
  void logmesg(CAC *, const std::string &);
  void expand_args(const char *, std::vector<std::string> &);
  void get_particle_pos(int, int, double &, double &, double &);

};

}

#endif

