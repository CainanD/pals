#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <limits.h>
#include <ncurses.h>

#include "db_parse.h"
#include "map.h"
#include "heap.h"

#define KEY_ESC 27
//#define KEY_UP 3
//#define KEY_DOWN 2
//#define KEY_ENTER 10

#define ROWS 21
#define COLS 80
#define WORLDSIZE 401

enum settings {
    NUMTRAINERS,
    SETTINGS_COUNT
};

enum directions {
    NORTH,
    SOUTH,
    EAST,
    WEST,
    NUM_DIRECTIONS
};

int settings[SETTINGS_COUNT] = {
    [NUMTRAINERS] = 6
};

void displayMap(int curX, int curY, map_t* world[WORLDSIZE][WORLDSIZE])
{
    int curRow = WORLDSIZE/2 - curY;
    int curCol = WORLDSIZE/2 + curX;

    if(!world[curRow][curCol])
    {
        int n;
        if(curRow <= 0 || (world[curRow-1][curCol]==NULL)){
            n = rand()%(COLS-2)+1;}
        else{
            n = world[curRow-1][curCol]->s;
        }

        int s;
        if(curRow >= WORLDSIZE-1 || (world[curRow+1][curCol]==NULL)){
            s = rand()%(COLS-2)+1;}
        else{
            s = world[curRow+1][curCol]->n;
        }

        int e; //Needs to be even for correct maze generation
        //printf("Value: %i", world[curRow][curCol+1]);
        if(curCol >= WORLDSIZE-1 || (world[curRow][curCol+1]==NULL)){
            e = 2*(rand()%((ROWS-2)/2))+1; 
        }
        else{
            //printMap(world[curRow][curCol+1]);
            e = world[curRow][curCol+1]->w;
        }

        int w;
        if(curCol <= 0 || (world[curRow][curCol-1]==NULL)){
            w = 2*(rand()%((ROWS-2)/2))+1;
        }
        else{
            w = world[curRow][curCol-1]->e;
        }

        int distance = abs(curX) + abs(curY) + 1;

        world[curRow][curCol] = generateMap(n, s, e, w, distance, settings[NUMTRAINERS]);

        //Edge cases
        checkEdges(curRow, curCol, world[curRow][curCol]);
    }
    printMap(world[curRow][curCol]);
}

map_t* getMap(int curX, int curY, map_t* world[WORLDSIZE][WORLDSIZE])
{
    return world[WORLDSIZE/2 - curY][WORLDSIZE/2 + curX];
}

void getSettings(int argc, char *argv[])
{
    for(int i = 0; i < argc; i++)
    {
        if(strcmp(argv[i], "--numtrainers") == 0)
        {
            settings[NUMTRAINERS] = atoi(argv[i+1]);
        }
    }
}

int* getFlyingCoords(map_t* map)
{
    echo();
    curs_set(TRUE);

    int* coords = (int*) malloc(2*sizeof(int)); 
    coords[0] = -1000;
    coords[1] = -1000;

    while(abs(coords[0]) > WORLDSIZE/2 || abs(coords[1]) > WORLDSIZE/2)
    {
        erase();
        printMap(map);
        printw("Fly Mode Activated! Enter two integers separated by a space.\n");
        printw("Enter the co-ordinates you want to fly to: ");
        scanw("%d %d", &coords[0], &coords[1]);
    }

    curs_set(0);
    noecho();

    return coords;
}

pokemon_t* welcome_player()
{
    pokemon_t* pal1 = createPokemon(1);
    pokemon_t* pal2 = createPokemon(1);
    pokemon_t* pal3 = createPokemon(1);
    pal1->level = 1;
    pal2->level = 1;
    pal3->level = 1;

    int index = 0;  
    int input = 0;

    while(input != 10)
    {
        if(input == KEY_UP && index > 0)   {index--;}
        if(input == KEY_DOWN && index < 2) {index++;}

        erase();
        printw("Welcome to A Mazing World of Pals! \n");
        printw("To start your journey, select one of the three pals using the ENTER button. \n");

        printInfo(pal1);
        printInfo(pal2);
        printInfo(pal3);

        mvprintw(4 + 4*index, 52, "<----");
        
        input = getch();
        refresh();
    }
    
    if(index == 0)
    {
        delete pal2;
        delete pal3;
        return pal1;
    }
    else if(index == 1)
    {
        delete pal1;
        delete pal3;
        return pal2;
    }
    else
    {
        delete pal1;
        delete pal2;
        return pal3;
    }
}

void init_terminal(void)
{
    initscr();
    raw();                 //unbuffered input
    noecho();              //does not print typed characters
    curs_set(0);           //Turns off cursor
    keypad(stdscr, TRUE);  //
    start_color();         //Makes terminal colorful
    initColors();
}

void gameOver()
{
    erase();
    attron(COLOR_PAIR(COLOR_RED));
    mvprintw(11, 35, "Game Over");
    attroff(COLOR_PAIR(COLOR_RED));
    refresh();
    getch();
}

/*
trainer_t* spawnPC(map_t* map)
{
    map->terrain[map->trainers[0].rowPos][map->trainers[0].colPos].occupied = FALSE;
    trainer_t pc_data = createTrainer(map, 0);
    trainer_t* pc = &(pc_data);
    return pc;
} */

int main(int argc, char *argv[])
{
    srand(time(NULL));
    init_terminal();
    db_parse(false);
    getSettings(argc, argv);
    
    map_t* world[WORLDSIZE][WORLDSIZE];
    for(int i = 0; i < WORLDSIZE; i++)
    {
        for(int j = 0; j < WORLDSIZE; j++)
        {
            world[i][j] = NULL;
        }
    }

    pokemon_t* starter_pokemon;
    starter_pokemon = welcome_player();

    int dir = 0;
    int curX = 0; 
    int curY = 0;
    int input;

    displayMap(curX, curY, world);
    map_t* map = getMap(curX, curY, world);
    //map->terrain[map->trainers[0].rowPos][map->trainers[0].colPos].occupied = FALSE;
    trainer_t pc_data = createTrainer(map, 0);
    trainer_t* pc = &(pc_data);
    movePC(map, dir, pc);
    pc->pokemon_team[0] = starter_pokemon;
    pc->numPokemon = 1;

    gameLoop:
    while(1)
    {
        displayMap(curX, curY, world);
        map = getMap(curX, curY, world);
        movePC(map, dir, pc);
        pc = &(map->trainers[0]);
        displayMap(curX, curY, world);
        trainer_t* trainers = map->trainers;
        MinHeap* q = map->q;

        while(1)
        {
            HeapNode* trainerNode = dequeue(q);
            trainer_t* trainer = (trainer_t*) trainerNode->o;
            moveTrainer(map, trainer);
            map->time = trainerNode->distance; //distance represents finish time

            if(trainer->type == 0) //if the trainer is the PC
            {
                if(map->terrain[trainer->rowPos][trainer->colPos].type == GRASS)
                {
                    if(rand()%10 == 0)
                    {
                        pokemon_t* otherPokemon = createPokemon(map->distance);
                        battlePokemon(pc, otherPokemon);
                    }
                }

                if(pc->defeated == TRUE) {goto Quit;}

                int i = 0;
                while(validNextMove(map, trainer) == FALSE || i == 0)
                {
                    erase();
                    displayMap(curX, curY, world);
                    printw("Coords: %d %d, PC location: %d %d, Current Time: %d\n", curX, curY, pc->rowPos, pc->colPos, map->time);
                    printw("Input {h,j,k,l} to move your character, f to fly, or Q to quit.\n");
                    refresh();

                    int nextRow = trainer->rowPos;
                    int nextCol = trainer->colPos;
                    input = getch();

                    switch(input)
                    {
                        case '1':
                        case 'b':
                        case 'z':
                            nextRow++;
                            nextCol--;
                            break;

                        case '2':
                        case 'j':
                        case 'x':
                        case KEY_DOWN:
                            nextRow++;
                            break;

                        case '3':
                        case 'n':
                        case 'c':
                            nextRow++;
                            nextCol++;
                            break;

                        case '4':
                        case 'h':
                        case 'a':
                        case KEY_LEFT:
                            nextCol--;
                            break;

                        case '5':
                        case ' ':
                        case '.':
                        case 's':
                            break;

                        case '6':
                        case 'l':
                        case 'd':
                        case KEY_RIGHT:
                            nextCol++;
                            break;

                        case '7':
                        case 'y':
                        case 'q':
                            nextRow--;
                            nextCol--;
                            break;

                        case '8':
                        case 'k':
                        case 'w':
                        case KEY_UP:
                    
                            nextRow--;
                            break;

                        case '9':
                        case 'u':
                        case 'e':
                            nextRow--;
                            nextCol++;
                            break;
                        
                        case '>':
                            if(map->terrain[nextRow][nextCol].symbol == 'C')
                            {
                                for(int i = 0; i < pc->numPokemon; i++)
                                {
                                    pc->pokemon_team[i]->hp = pc->pokemon_team[i]->stats[HP];
                                }

                                erase();
                                printw("All of your pokemon have been healed. \n");
                                printw("Have a fantastic day! \n");
                                printw("Input < to exit. \n");
                                refresh();

                                while(input != '<' && input != KEY_ESC)
                                {
                                    input = getch();
                                }
                            }
                            else if(map->terrain[nextRow][nextCol].symbol == 'M')
                            {
                                pc->bag[POTIONS] = 3;
                                pc->bag[REVIVES] = 2;
                                pc->bag[POKEBALLS] = 5;

                                erase();
                                printw("All of your items have been restored. \n");
                                printw("Have a fantastic day! \n");
                                printw("Input < to exit. \n");
                                refresh();

                                while(input != '<' && input != KEY_ESC)
                                {
                                    input = getch();
                                }
                            }

                            i = -1;
                            break;
                        
                        case 't':
                            {
                                input = ' ';
                                int startIndex = 0;
                                int scrollDir = 0;

                                while(input != 't' && input != KEY_ESC)
                                {
                                    if(scrollDir == KEY_UP && startIndex > 0) {startIndex--;}
                                    else if(scrollDir == KEY_DOWN && 
                                    startIndex < settings[NUMTRAINERS] - ROWS) {startIndex++;}

                                    printTrainers(map, startIndex);
                                    scrollDir = getch();
                                    input = scrollDir;
                                }
                            }
                            i = -1;
                            break;

                        case 'B':
                        {
                            erase();
                            int bagIndex = enterBag(pc);

                            if(bagIndex == POTIONS)
                            {
                                if(!healPokemon(pc))
                                {
                                    pc->bag[POTIONS]++;
                                }
                            }
                            else if(bagIndex == REVIVES)
                            {
                                if(revivePokemon(pc) == -1)
                                {
                                    pc->bag[REVIVES]++;
                                }
                            }
                            else if(bagIndex == POKEBALLS)
                            {
                                pc->bag[POKEBALLS]++;
                            }
                            break;
                        }

                        case 'P':
                            erase();
                            for(int i = 0; i < pc->numPokemon; i++)
                            {   
                                printInfo(pc->pokemon_team[i]);
                            }

                            refresh();
                            getch();
                            break;


                        case 'Q':
                        case 'm':
                            goto Quit;
                            break;

                        case 'f':
                        {   
                            int* coords = getFlyingCoords(map);
                            curX = coords[0];
                            curY = coords[1];
                        }
                            trainerNode->distance = map->time;
                            enqueue(q, trainerNode);

                            goto gameLoop;
                            break;


                        default:
                            i = -1;
                            break;
                    }

                    trainer->nextCol = nextCol;
                    trainer->nextRow = nextRow;
                    i++;
                }

                #define GATE 3
                if(validNextMove(map, trainer) == GATE) //If next move is gate
                {
                    if(trainer->nextRow == 0) {dir = NORTH; curY++;}
                    else if(trainer->nextRow == ROWS - 1) {dir = SOUTH; curY--;}
                    else if(trainer->nextCol == COLS - 1) {dir = EAST; curX++;}
                    else if(trainer->nextCol == 0) {dir = WEST; curX--;}
                    trainerNode->distance = map->time;
                    enqueue(q, trainerNode);
                    break;
                }

                //Sorry about the mouthful
                trainerNode->distance = map->time + trainer->movement[map->terrain[trainer->nextRow][trainer->nextCol].type];
                enqueue(q, trainerNode);
            }
            else
            {
                int** distanceMap = calcDistanceFromPC(map, trainer, &trainers[0]);
                int finishTime = calcNextMove(map, trainer, distanceMap, map->time);
                freeDistanceMap(distanceMap);
                trainerNode->distance = finishTime;
                enqueue(q, trainerNode);
            }
        }
    }

    Quit:
    gameOver();

    /*
    for(int i = 0; i < WORLDSIZE; i++)
    {
        for(int j = 0; j < WORLDSIZE; j++)
        {
            if(world[i][j])
            {
                freeMap(world[i][j]);
            }
        }
    } I like it leaky */

    endwin();
    return 0;
}