#ifndef MAP_H
#define MAP_H

#include "heap.h"
#include "pokemon.h"
#define TEAMSIZE 6

typedef struct terrain
{
    int type;
    int occupied; 
    char symbol;
    short color;

} terrain_t;

enum terrainTypes {
    WALL,     // '%'
    PASSAGE,  // ':' temporary for maze building
    DEADEND,  // ':' important for water generation and pathfinding
    GRASS,    // ':'
    ROAD,     // '#'
    ENTRANCE, // '#' special case of road
    WATER,    // '~'
    MART,     // 'M'
    CENTER,   // 'C'
    CLEARING, // '.'
    TREE,     // 'T'
    NUM_TERRAIN_TYPES
};

enum bagItems {
    POTIONS,
    REVIVES,
    POKEBALLS,
    NUM_ITEM_TYPES
};

typedef class trainer
{
    public:
    int type;
    int rowPos;
    int colPos;
    char symbol;
    int* movement;
    int defeated;
    short color;

    int numPokemon;
    pokemon_t* pokemon_team[TEAMSIZE];
    int* bag;

    int nextRow;
    int nextCol;

    int direction;
    terrain_t spawnTerrain;
} trainer_t;

enum trainerTypes {
    PC,
    RIVAL,
    HIKER,
    PACER,
    WANDERER,
    SENTRY,
    EXPLORER,
    NUM_TRAINER_TYPES
};

typedef class map
{
    public:
    int n;
    int s;
    int e;
    int w;
    int rows;
    int cols;
    int numTrainers;
    int distance;
    terrain_t** terrain;
    trainer_t* trainers;

    int time;
    MinHeap* q;
    HeapNode* trainerNodes;
} map_t;

void initColors(void);
void printMap(map_t* map);
map_t* generateMap(int n, int s, int e, int w, int distance, int numTrainers);
void freeMap(map_t* map);
void checkEdges(int curRow, int curCol, map_t* map);

int** calcDistanceFromPC(map_t* map, trainer_t* trainer, trainer_t* pc);
int calcNextMove(map_t* map, trainer_t* trainer, int** distanceMap, int time);
void moveTrainer(map_t* map, trainer_t* trainer);
int validNextMove(map_t* map, trainer_t* trainer);
void freeDistanceMap(int** distanceMap);

void battlePokemon(trainer_t* pc, pokemon_t* pokemon);
pokemon_t* swapPokemon(trainer_t* pc);
int revivePokemon(trainer_t* pc);
int healPokemon(trainer_t* pc);
int enterBag(trainer_t* pc);
pokemon_t* getAlivePokemon(trainer_t* trainer);
void battlePokemon(trainer_t* pc, pokemon_t* pokemon);


void movePC(map_t* map, int dir, trainer_t* pc);
trainer_t createTrainer(map_t* map, int type);
void spawnTrainers(map_t* map, int numTrainers);
void printTrainers(map_t* map, int startIndex);
void destroyTrainer(trainer_t* trainer);

#endif