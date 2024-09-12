#ifndef PAL_H
#define PAL_H
#include <string>
#include <vector>
#include "db_parse.h"

enum IVs {
    HP,
    ATTACK,
    DEFENSE,
    SPEED,
    SPECIAL_ATTACK,
    SPECIAL_DEFENSE,
    NUM_IVS
};

typedef class pokemon
{
    public:
    int id;
    char name[30];
    int level;
    int moves[2];
    bool isShiny;
    int type;
    int hp;
    int baseSpeed;
    
    int ivs[NUM_IVS];
    int stats[NUM_IVS];

    pokemon() : id(-1) {}
} pokemon_t;

pokemon_t* createPokemon(int distance);
void destroyPokemon(pokemon_t* pokemon);
void printInfo(pokemon_t* pokemon);

std::string attack(pokemon_t* agressor, int moveIndex, pokemon_t* defender);

#endif