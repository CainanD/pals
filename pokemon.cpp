#include <stdlib.h>
#include <ncurses.h>
#include <cmath>
#include <cstring>
#include <sstream>
#include <limits.h>

#include "pokemon.h"
#include "map.h"

#define STRUGGLE 165
#define KEY_ESC 27


void printInfo(pokemon_t* pal)
{
    if(pal->hp <= 0)
    {
        attron(COLOR_PAIR(8));
    }
    else if (pal->isShiny)
    {
        attron(COLOR_PAIR(COLOR_CYAN));
    }

    printw("\nLevel %d %s, %d/%d HP\n", pal->level, pal->name, pal->hp, pal->stats[HP]);
    printw("Stats: %d Atk, %d Def, %d SpA, %d SpD, %d Spe \n", pal->stats[ATTACK],
    pal->stats[DEFENSE], pal->stats[SPECIAL_ATTACK], pal->stats[SPECIAL_DEFENSE], pal->stats[SPEED]);
    printw("Moves: %s %s \n", moves[pal->moves[0]].identifier, moves[pal->moves[1]].identifier);

    
    if(pal->hp <= 0)
    {
        attroff(COLOR_PAIR(8));
    }
    else if (pal->isShiny)
    {
        attroff(COLOR_PAIR(COLOR_CYAN));
    }
}

pokemon_t* createPokemon(int distance)
{
    pokemon_t* pal = new pokemon_t;
    int id = pokemon[rand()%898].id;
    int level;

    if(distance > 200)
    {
        level = rand()%100 + 1;
    }
    else
    {
        level = rand()%((distance+1)/2) + 1;
    }


    std::vector<int> possibleMoves;
    int i = 0;

    while(pokemon_moves[i].pokemon_id != id) {i++;}

    while(possibleMoves.size() == 0)
    {
        int j = 0;
        while(pokemon_moves[i+j].pokemon_id == id)
        {
            if(pokemon_moves[i+j].level <= level && pokemon_moves[i+j].pokemon_move_method_id == 1)
            {
                possibleMoves.push_back(pokemon_moves[i+j].move_id);
            }
            j++;
        }
        level++;
    }

    if(possibleMoves.size() == 1)
    {
        pal->moves[0] = possibleMoves[0];
        pal->moves[1] = STRUGGLE;
    }
    else
    {
        int index = rand()%possibleMoves.size();
        pal->moves[0] = possibleMoves[index];
        int newIndex = rand()%possibleMoves.size();
        while(newIndex == index) {newIndex = rand()%possibleMoves.size();}
        pal->moves[1] = possibleMoves[newIndex];
    }

    int type;
    for(int i = 0; i < 1673; i++)
    {
        if(pokemon_types[i].pokemon_id == id && pokemon_types[i].slot == 1)
        {
            type = pokemon_types[i].type_id;
        }
    }

    for(int k = 0; k < NUM_IVS; k++)
    {
        pal->ivs[k] = rand()%16;

        if(k == HP)
        {
            pal->stats[HP] = (2*level*(pokemon_stats[NUM_IVS*(id - 1) + k + 1].base_stat + pal->ivs[HP])/100) + level + 10;
        }
        else
        {
            pal->stats[k] = (2*level*(pokemon_stats[NUM_IVS*(id - 1) + k + 1].base_stat + pal->ivs[k])/100) + 5;

            if (k == SPEED)
            {
                pal->baseSpeed = pokemon_stats[NUM_IVS*(id - 1) + k + 1].base_stat;
            }
        }
    }

    bool isShiny = (rand()%10 == 0);

    pal->type = type;
    pal->hp = pal->stats[HP];
    pal->id = id;
    strcpy(pal->name, pokemon[id].identifier);
    pal->isShiny = isShiny;
    pal->level = level;

    return pal;
}

std::string attack(pokemon_t* agressor, int moveIndex, pokemon_t* defender)
{
    std::stringstream ss;
    move_db move = moves[agressor->moves[moveIndex]];

    ss << agressor->name << " uses " << move.identifier << "." << std::endl;

    if(rand()%100 > move.accuracy)
    {
        ss << "The move misses..." << std::endl;
        return ss.str();
    }

    float stab =  1;
    float random = float(100 - rand()%15)/100.0;
    float critical = 1;

    if(rand()%256 < agressor->baseSpeed/2)
    {
        critical = 1.5;
        ss << "A critical hit!" << std::endl;
    }

    if(agressor->type == move.type_id)
    {
        stab = 1.5;
        ss << "The move is super effective!" << std::endl;
    }

    int damageDone;
    if(move.power == INT_MAX)
    {
        damageDone = 1;
    }
    else
    {
        float damage = (((2/5)*agressor->level+2)*move.power*(agressor->stats[ATTACK]/defender->stats[DEFENSE])/50 + 2) * critical * stab * random;
        damageDone = std::max(1, (int)damage);
    }
    defender->hp -= damageDone;
    ss << "The move does " << damageDone << " damage." << std::endl;

    return ss.str();
}

void destroyPokemon(pokemon_t* pokemon)
{
    delete pokemon;
}