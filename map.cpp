#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>
#include <ncurses.h>
#include <string.h>
#include <sstream>

#include "map.h"
#include "heap.h"

#define FALSE 0
#define TRUE  1
#define KEY_ESC 27

#define WORLDSIZE 401
#define ROWS 21
#define COLS 80

const char trainerSymbols[NUM_TRAINER_TYPES] = {
    '@',
    'r',
    'h',
    'p',
    'w',
    's',
    'e'
};

const char terrainSymbols[NUM_TERRAIN_TYPES] = {
    '%',
    ':',
    ':',
    ':',
    '#',
    '#',
    '~',
    'M',
    'C',
    '.',
    'T'
};

void initColors(void)
{
    #define COLOR_GRAY 8
    init_pair(COLOR_WHITE, COLOR_WHITE, COLOR_BLACK);
    init_pair(COLOR_RED,  COLOR_RED,  COLOR_BLACK);
    init_pair(COLOR_GREEN, COLOR_GREEN, COLOR_BLACK);
    init_pair(COLOR_YELLOW, COLOR_YELLOW, COLOR_BLACK);
    init_pair(COLOR_BLUE, COLOR_BLUE, COLOR_BLACK);
    init_pair(COLOR_CYAN, COLOR_CYAN, COLOR_BLACK);
    init_pair(COLOR_GRAY, COLOR_GRAY, COLOR_BLACK);
}

const short terrainColors[NUM_TERRAIN_TYPES] = {
    COLOR_PAIR(COLOR_GRAY),
    COLOR_PAIR(COLOR_GREEN),
    COLOR_PAIR(COLOR_GREEN),
    COLOR_PAIR(COLOR_GREEN),
    COLOR_PAIR(COLOR_YELLOW),
    COLOR_PAIR(COLOR_YELLOW),
    COLOR_PAIR(COLOR_BLUE),
    COLOR_PAIR(COLOR_CYAN),
    COLOR_PAIR(COLOR_RED),
    COLOR_PAIR(COLOR_GREEN),
    COLOR_PAIR(COLOR_GREEN)
};

enum directions {
    NORTH,
    SOUTH,
    EAST,
    WEST,
    NUM_DIRECTIONS
};

void add(int row, int col, int arr[], int cols)
{
    int end = arr[0] + 1;
    int element = cols*row+col;
    int containsElement = 0;
    for(int i = 1; i <= end; i++)
    {
        if(arr[i] == element)
        {
            containsElement = 1;
        }
    }

    if(containsElement == 0)
    {
        arr[end] = element;
        arr[0]++;
    }
}

int contains(int arr[], int k)
{
    int i = 1;
    while(arr[i] != -1 && i < ROWS*COLS)
    {
        if (arr[i] == k)
        {
            return 1;
        }
        i++;
    }
    return 0;
}

int pickRandom(int* arr)
{
    int index = rand()%arr[0] + 1;
    int picked = arr[index];

    int i = index;
    while(arr[i] != -1) 
    {
        arr[i] = arr[i+1];
        i++;
    }

    arr[0]--;
    return picked;
}

int min(int a, int b)
{
    if(a < b)
    {
        return a;
    }
    else
    {
        return b;
    }
}

int** generateMaze(int n, int s, int e, int w, int rows, int cols)
{
    int** maze = (int**)malloc(rows*sizeof(int*));
    for(int i = 0; i < rows; i++)
    {
        maze[i] = (int*)malloc(cols*sizeof(int));
        for(int j = 0; j < cols; j++)
        {
            maze[i][j] = WALL;
        }
    }

    //Prim's Algorithm (from Wikipedia)
    int frontierList[rows*cols];
    for(int i = 0; i < rows*cols; i++)
    {
        frontierList[i] = -1;
    }

    //Element at first index indicates size of list
    frontierList[0] = 0;

    //Initializing maze at entrance
    add(w-1, 0, frontierList, cols);
    maze[w-1][0] = PASSAGE;

    while(frontierList[0] > 0)
    {
        int index = pickRandom(frontierList);
        int wallRow = index / cols;
        int wallCol = index % cols;
        maze[wallRow][wallCol] = PASSAGE;

        int passages[6] = {0, -1, -1, -1, -1, -1};

        if(wallRow >= 2)
        {
            if (maze[wallRow-2][wallCol] == WALL)
            {
                add(wallRow-2, wallCol, frontierList, cols);
            }
            else
            {
                add(wallRow-1, wallCol, passages, cols);
            }
        }
        if(wallRow < rows - 2)
        {
            if (maze[wallRow+2][wallCol] == WALL)
            {
                add(wallRow+2, wallCol, frontierList, cols);
            }
            else
            {
                add(wallRow+1, wallCol, passages, cols);
            }
        }
        if(wallCol >= 2)
        {
            if (maze[wallRow][wallCol-2] == WALL)
            {
                add(wallRow, wallCol-2, frontierList, cols);
            }
            else
            {
                add(wallRow, wallCol-1, passages, cols);
            }
        }
        if(wallCol < cols - 2)
        {
            if (maze[wallRow][wallCol+2] == WALL)
            {
                add(wallRow, wallCol+2, frontierList, cols);
            }
            else
            {
                add(wallRow, wallCol+1, passages, cols);
            }
        }

        if (passages[0] > 0)
        {
            int passageIndex = pickRandom(passages);
            int passageRow = passageIndex / cols;
            int passageCol = passageIndex % cols;
            maze[passageRow][passageCol] = PASSAGE;
        }
    }

    //Opening up passages to the gates
    maze[e-1][cols-1] = PASSAGE;
    maze[w-1][0] = PASSAGE;
    maze[e-1][cols-2] = PASSAGE;
    maze[0][n-1] = PASSAGE;
    maze[rows-1][s-1] = PASSAGE;

    return maze;
}

void solveMazeRec(int dest, int** maze, int rows, int cols, int* path)
{
    int pathLength = path[0];
    if (path[pathLength] == dest)
    {
        int startIndex = path[1];
        int startRow = startIndex / cols;
        int startCol = startIndex % cols;

        int destRow = dest / cols;
        int destCol = dest % cols;

        //Prevents duplicate paths
        if(maze[startRow][startCol] == ENTRANCE 
        && maze[destRow][destCol] == ENTRANCE)
        {
            return;
        }

        maze[startRow][startCol] = ENTRANCE;
        for (int i = 2; i < pathLength; i++)
        {
            int curRow = path[i] / cols;
            int curCol = path[i] % cols;
            if(maze[curRow][curCol] != ENTRANCE)
            {
                maze[curRow][curCol] = ROAD;
            }
        }
        maze[destRow][destCol] = ENTRANCE;
        
        return;
    }
    else
    {
        int curIndex = path[path[0]];
        int curCol = curIndex % cols;
        int curRow = curIndex / cols;
        int rowMoves[4] = {-1, 1, 0, 0};
        int colMoves[4] = {0, 0, -1, 1};

        int validMoves = 0;
        for(int i = 0; i < 4; i++)
        {
            int newRow = curRow + rowMoves[i];
            int newCol = curCol + colMoves[i];
            int newIndex = (newRow)*cols + newCol;

            if(newRow >= 0 && newRow < rows && newCol >= 0 && newCol < cols 
            && contains(path, newIndex) == 0 && maze[newRow][newCol] != WALL)
            {
                validMoves++;
                int newPath[rows*cols];
                for(int i = 0; i < rows*cols; i++)
                {
                    newPath[i] = path[i];
                }
                newPath[0]++;
                newPath[newPath[0]] = newIndex;
                solveMazeRec(dest, maze, rows, cols, newPath);
            }
        }
        if(validMoves == 0)
        {
            if(maze[curRow][curCol] != ROAD && maze[curRow][curCol] != ENTRANCE)
            {
                maze[curRow][curCol] = DEADEND; //Used for generating water
            }
            return;
        }
    }
}

void solveMaze(int src, int dest, int rows, int cols, int** maze)
{
    int path[ROWS*COLS];
    for(int i = 0; i < ROWS*COLS; i++)
    {
        path[i] = -1;
    }

    path[0] = 1;
    path[1] = src;
    solveMazeRec(dest, maze, rows, cols, path);

    return; 
}

void placeWater(int rows, int cols, int** maze)
{
    int placedWater = 0;
    for(int i = 1; i < rows - 1; i++)
    {
        for(int j = 1; j < cols - 1; j++)
        {
            int wallCount = 0;
            for(int iOffset = -1; iOffset <= 1; iOffset++)
            {
                for(int jOffset = -1; jOffset <= 1; jOffset++)
                {
                    if(maze[i+iOffset][j+jOffset] == WALL 
                    || maze[i+iOffset][j+jOffset] == DEADEND)
                    {
                        wallCount++;
                    }
                }
            }

            if(wallCount == 9)
            {
                placedWater++;
                for(int iOffset = -1; iOffset <= 1; iOffset++)
                {
                    for(int jOffset = -1; jOffset <= 1; jOffset++)
                    {
                        maze[i+iOffset][j+jOffset] = WATER;
                    }
                }
            }
        }
    }
}

void placeClearings(int rows, int cols, int** maze, int distance)
{
    int clearingsPlaced = 0;
    int centerPlaced = 0;
    int martPlaced = 0;

    //Determines the distribution of centers and marts
    //as a function of distance from spawn
    if(distance > 200)
    {
        if(rand()%20 != 0)
        {
            martPlaced = 1;
        }
        if(rand()%20 != 0)
        {
            centerPlaced = 1;
        }
    }
    else if (distance > 20)
    {
        //Buildings have 50% chance of spawning at a distance of 20
        //Buildings have a 20% chance of spawning at a distance of 100
        if(rand()%(distance/20) != 0) 
        {
            martPlaced = 1;
        }
        if(rand()%(distance/20) != 0) 
        {
            centerPlaced = 1;
        }
    }

    for(int i = 0; i < rows; i++)
    {
        for(int j = 0; j < cols; j++)
        {
            int numRoads = 0;
            if(maze[i][j] == ROAD || maze[i][j] == ENTRANCE)
            {
                if(maze[i][j] == ENTRANCE)
                {
                    numRoads++;
                }
                if(i-1 > 0)
                {
                    if(maze[i-1][j] == ROAD || maze[i-1][j] == ENTRANCE){numRoads++;}
                }
                if(j-1 > 0)
                {
                    if(maze[i][j-1] == ROAD || maze[i][j-1] == ENTRANCE){numRoads++;}
                }
                if(i+1 < rows)
                {
                    if(maze[i+1][j] == ROAD || maze[i+1][j] == ENTRANCE){numRoads++;}
                }
                if(j+1 < cols)
                {
                    if(maze[i][j+1] == ROAD || maze[i][j+1] == ENTRANCE){numRoads++;}
                }
            }

            if(numRoads > 2)
            {
                clearingsPlaced++;
                int clearingSize = 3;
                for(int k = -clearingSize; k <= clearingSize; k++)
                {
                    for(int l = -clearingSize; l <= clearingSize; l++)
                    {
                        if(i+k >= 0 && i+k < rows && j+l >= 0 && j+l < cols
                        && maze[i+k][j+l] != ROAD && maze[i+k][j+l] != ENTRANCE
                        && maze[i+k][j+l] != MART && maze[i+k][j+l] != CENTER)
                        {
                            maze[i+k][j+l] = CLEARING;
                        }
                    }
                }

                //Placing buildings
                for(int k = -(clearingSize-1); k <= clearingSize; k++)
                {
                    for(int l = -(clearingSize-1); l <= clearingSize; l++)
                    {
                        if(centerPlaced == 0 && i+k >= 0 && i+k+1 < rows && j+l >= 0 && j+l+1 < cols
                        && maze[i+k][j+l] == CLEARING && maze[i+k][j+l+1] == CLEARING
                        && maze[i+k+1][j+l] == CLEARING && maze[i+k+1][j+l+1] == CLEARING)
                        {
                            centerPlaced = 1;
                            maze[i+k][j+l] = CENTER;
                            maze[i+k+1][j+l] = CENTER;
                            maze[i+k][j+l+1] = CENTER;
                            maze[i+k+1][j+l+1] = CENTER;
                        }

                        if(martPlaced == 0 && i-k-1 >= 0 && i-k < rows && j-l-1 >= 0 && j-l < cols
                        && maze[i-k][j-l] == CLEARING && maze[i-k-1][j-l] == CLEARING
                        && maze[i-k][j-l-1] == CLEARING && maze[i-k-1][j-l-1] == CLEARING)
                        {
                            martPlaced = 1;
                            maze[i-k][j-l] = MART;
                            maze[i-k-1][j-l] = MART;
                            maze[i-k][j-l-1] = MART;
                            maze[i-k-1][j-l-1] = MART;
                        }
                    }
                }
            }
        }
    }
}

void freeMaze(int** maze, int rows)
{
    for(int i = 0; i < rows; i++)
    {
        free(maze[i]);
    }
    free(maze);
}

terrain_t createTerrain(int type)
{
    terrain_t terrain;
    terrain.occupied = 0;
    terrain.type = type;
    terrain.symbol = terrainSymbols[type];
    terrain.color = terrainColors[type];
    
    return terrain;
}

map_t* generateMap(int n, int s, int e, int w, int distance, int numTrainers)
{
    //Initialize the map of completely boulders
    terrain_t** terrain = (terrain_t**)malloc(ROWS*sizeof(terrain_t*));
    for(int i = 0; i < ROWS; i++)
    {
        terrain[i] = (terrain_t*)malloc(COLS*sizeof(terrain_t));
        for(int j = 0; j < COLS; j++)
        {
            terrain[i][j] = createTerrain(WALL);
        }
    }

    terrain[0][n] =      createTerrain(ENTRANCE);
    terrain[ROWS-1][s] = createTerrain(ENTRANCE);
    terrain[w][0] =      createTerrain(ENTRANCE);
    terrain[e][COLS-1] = createTerrain(ENTRANCE);

    //Create the maze
    int mazeRows = ROWS - 2;
    int mazeCols = COLS - 2;
    int** maze = generateMaze(n, s, e, w, mazeRows, mazeCols);

    //Solve the maze for roads
    int srcN = n-1;
    int destS = (mazeCols)*(mazeRows-1)+(s-1);
    int srcW = (w-1)*(mazeCols);
    int destE = (e-1)*(mazeCols) + (mazeCols - 1);
    solveMaze(srcN, destS, mazeRows, mazeCols, maze);
    solveMaze(srcW, destE, mazeRows, mazeCols, maze);

    //Placing water
    placeWater(mazeRows, mazeCols, maze);

    //Placing clearings
    placeClearings(mazeRows, mazeCols, maze, distance);

    for (int i = 0; i < mazeRows; i++)
    {
        for (int j = 0; j < mazeCols; j++)
        {
            if(maze[i][j] == WALL)
            {
                terrain[i+1][j+1] = createTerrain(TREE);
            }
            if(maze[i][j] == PASSAGE)
            {
                terrain[i+1][j+1] = createTerrain(GRASS);
            }
            if(maze[i][j] == DEADEND)
            {
                terrain[i+1][j+1] = createTerrain(DEADEND);
            }
            if(maze[i][j] == ROAD || maze[i][j] == ENTRANCE)
            {
                terrain[i+1][j+1] = createTerrain(ROAD);
            }
            if(maze[i][j] == WATER)
            {
                terrain[i+1][j+1] = createTerrain(WATER);
            }
            if(maze[i][j] == CLEARING)
            {
                terrain[i+1][j+1] = createTerrain(CLEARING);
            }
            if(maze[i][j] == CENTER)
            {
                terrain[i+1][j+1] = createTerrain(CENTER);
            }
            if(maze[i][j] == MART)
            {
                terrain[i+1][j+1] = createTerrain(MART);
            }
        }
    }

    freeMaze(maze, mazeRows);
    
    map_t* map = (map_t*) malloc(sizeof(map_t));
    map->n = n;
    map->s = s;
    map->e = e;
    map->w = w;
    map->rows = ROWS;
    map->cols = COLS;
    map->terrain = terrain;
    map->time = 0;
    map->distance = distance;

    spawnTrainers(map, numTrainers);
    trainer_t* trainers = map->trainers;
    trainer_t* pc = &(trainers[0]);

    MinHeap* q = createHeap(ROWS*COLS);
    HeapNode* trainerNodes = (HeapNode*)malloc((numTrainers+1)*sizeof(HeapNode)); 

    for(int i = 0; i <= numTrainers; i++)
    {   
        int finishTime;
        trainer_t trainer = trainers[i];
        if(i == 0)
        {
            int nextRow = pc->rowPos;
            int nextCol = pc->colPos;
            trainer.nextRow = nextRow;
            trainer.nextCol = nextCol;
            finishTime = map->time;
        }
        else
        {
            int** distanceMap = calcDistanceFromPC(map, &trainer, pc);
            finishTime = calcNextMove(map, &trainer, distanceMap, map->time);
            freeDistanceMap(distanceMap);
        }
        
        trainerNodes[i].o = (void*) &trainers[i];
        trainerNodes[i].distance = finishTime;
        enqueue(q, &trainerNodes[i]);
    }

    map->q = q;
    map->trainerNodes = trainerNodes;

    return map;
}

void freeMap(map_t* map)
{
    for(int i = 0; i < ROWS; i++)
    {
        free(map->terrain[i]);
    }

    for(int i = 0; i < map->numTrainers; i++)
    {
        destroyTrainer(&(map->trainers[i]));
    }

    freeHeap(map->q);
    free(map->trainers);
    free(map->trainerNodes);
    free(map->terrain);
    free(map);
}

void checkEdges(int curRow, int curCol, map_t* map)
{
    if(curRow == 0)
    {
        for(int i = 0; i < COLS; i++)
        {
            map->terrain[0][i] = createTerrain(WALL);
        }
    }

    if(curRow == WORLDSIZE-1)
    {
        for(int i = 0; i < COLS; i++)
        {
            map->terrain[ROWS-1][i] = createTerrain(WALL);
        }
    }

    if(curCol == 0)
    {
        for(int i = 0; i < ROWS; i++)
        {
            map->terrain[i][0] = createTerrain(WALL);
        }
    }

    if(curCol == WORLDSIZE-1)
    {
        for(int i = 0; i < ROWS; i++)
        {
            map->terrain[i][COLS-1] = createTerrain(WALL);
        }
    }
}

int pcMovement[NUM_TERRAIN_TYPES] = {
    /*[WALL]     =*/ INT_MAX,
    /*[PASSAGE]  =*/ 20,
    /*[DEADEND]  =*/ 20,
    /*[GRASS]    =*/ 20,
    /*[ROAD]     =*/ 10,   
    /*[ENTRANCE] =*/ 10, 
    /*[WATER]    =*/ INT_MAX,
    /*[MART]     =*/ 10,
    /*[CENTER]   =*/ 10,
    /*[CLEARING] =*/ 10,                    
    /*[TREE]     =*/ INT_MAX
};

int defaultMovement[NUM_TERRAIN_TYPES] = {
    /*[WALL]     =*/ INT_MAX,
    /*[PASSAGE]  =*/ 20,
    /*[DEADEND]  =*/ 20,
    /*[GRASS]    =*/ 20,
    /*[ROAD]     =*/ 10,   
    /*[ENTRANCE] =*/ INT_MAX, 
    /*[WATER]    =*/ INT_MAX,
    /*[MART]     =*/ 50,
    /*[CENTER]   =*/ 50,
    /*[CLEARING] =*/ 10,                    
    /*[TREE]     =*/ INT_MAX
};

int hikerMovement[NUM_TERRAIN_TYPES] = {
    /*[WALL]     =*/ INT_MAX,
    /*[PASSAGE]  =*/ 15,
    /*[DEADEND]  =*/ 15,
    /*[GRASS]    =*/ 15,
    /*[ROAD]     =*/ 10,   
    /*[ENTRANCE] =*/ INT_MAX, 
    /*[WATER]    =*/ INT_MAX,
    /*[MART]     =*/ 50,
    /*[CENTER]   =*/ 50,
    /*[CLEARING] =*/ 10,                    
    /*[TREE]     =*/ 15
};

int* trainerMovement[NUM_TRAINER_TYPES] = {
    /*[PC]       =*/ pcMovement,
    /*[RIVAL]    =*/ defaultMovement,
    /*[HIKER]    =*/ hikerMovement,
    /*[PACER]    =*/ defaultMovement,
    /*[WANDERER] =*/ defaultMovement,
    /*[SENTRY]   =*/ defaultMovement,
    /*[EXPLORER] =*/ defaultMovement
}; 

/*
const int swimmerMovement[128] = {
    ['%'] = INT_MAX,
    ['T'] = INT_MAX,
    ['#'] = INT_MAX,
    ['C'] = INT_MAX,
    ['M'] = INT_MAX,
    [':'] = INT_MAX,
    ['.'] = INT_MAX,
    ['~'] = 7,
}; Unused NPC type */ 

int** calcDistanceFromPC(map_t* map, trainer_t* trainer, trainer_t* pc)
{
    int mapRows = map->rows;
    int mapCols = map->cols;
    int srcCol = pc->colPos;
    int srcRow = pc->rowPos;

    //Dual heap/map data structure allows updating distance values within the heap
    MinHeap* queue = createHeap(mapRows*mapCols);
    HeapNode** nodes = (HeapNode**)malloc(mapRows*sizeof(HeapNode*)); 
    
    for(int i = 0; i < mapRows; i++)
    {
        nodes[i] = (HeapNode*)malloc(mapCols*sizeof(HeapNode));
    }

    for(int i = 0; i < mapRows; i++)
    {
        for(int j = 0; j < mapCols; j++)
        {
            nodes[i][j].row = i;
            nodes[i][j].col = j;

            if(i == srcRow && j == srcCol)
            {
                nodes[i][j].distance = 0;
            }
            else
            {
                nodes[i][j].distance = INT_MAX;
            }

            enqueue(queue, &nodes[i][j]);
        }
    }

    while (queue->size != 0)
    {
        //printHeap(q);
        HeapNode* u = dequeue(queue);
        int row = u->row;
        int col = u->col;

        for(int i = -1; i <= 1; i++)
        {
            for(int j = -1; j <= 1; j++)
            {
                if(abs(i)+abs(j)!=0 && row+i>0 && row+i<mapRows-1 && col+j>0 && col+j<mapCols-1)
                {
                    HeapNode* v = &nodes[row+i][col+j];
                    int movementCost = trainer->movement[map->terrain[row+i][col+j].type];
                    //printf("Terrain: %c, Index: %d\n", map->terrain[row+i][col+j], map->terrain[row+i][col+j]);
                    //printf("Distance: %d, Movement Cost: %d\n", u->distance, movementCost);
                    if(u->distance != INT_MAX && movementCost != INT_MAX && u->distance + movementCost < v->distance)
                    {
                        modifyDistance(queue, v, u->distance + movementCost);
                    }
                }
            }
        }
    }

    int** distance = (int**)malloc(mapRows*sizeof(int*));
    for(int i = 0; i < mapRows; i++)
    {
        distance[i] = (int*)malloc(mapCols*sizeof(int));
        for(int j = 0; j < mapCols; j++)
        {
            distance[i][j] = nodes[i][j].distance;
            /**
            if(distance[i][j] != INT_MAX)
            {
                printf(" %02d", nodes[i][j].distance%100);
            }
            else
            {
                printf("   ");
            } */
        }
        //printf("\n");
    }

    for(int i = 0; i < mapRows; i++)
    {
        free(nodes[i]);
    }

    free(nodes);
    freeHeap(queue);

    return distance;
}

void freeDistanceMap(int** distanceMap)
{
    for(int i = 0; i < ROWS; i++)
    {
        free(distanceMap[i]);
    }

    free(distanceMap);
}

void movePC(map_t* map, int dir, trainer_t* pc_clone)
{
    map->terrain[map->trainers[0].rowPos][map->trainers[0].colPos].occupied = FALSE;
    map->trainers[0] = *pc_clone;
    trainer_t* pc = &map->trainers[0];

    if(dir == SOUTH)
    {
        pc->rowPos = 1;
        pc->colPos = map->n;
    }
    else if(dir == NORTH)
    {
        pc->rowPos = map->rows-2;
        pc->colPos = map->s;
    }
    else if(dir == WEST)
    {
        pc->colPos = map->cols-2;
        pc->rowPos = map->e;
    }
    else if(dir == EAST)
    {
        pc->colPos = 1;
        pc->rowPos = map->w;
    }

    if(pc->movement[map->terrain[pc->rowPos][pc->colPos].type] == INT_MAX
    || map->terrain[pc->rowPos][pc->colPos].occupied == TRUE)
    {
        int valid = FALSE;
        for(int i = -1; i <= 1; i++)
        {
            for(int j = -1; j <= 1; j++)
            {
                if(pc->movement[map->terrain[pc->rowPos][pc->colPos].type] != INT_MAX
                && map->terrain[pc->rowPos][pc->colPos].occupied == FALSE && valid == FALSE)
                {
                    pc->rowPos = pc->rowPos + i;
                    pc->colPos = pc->colPos + j;
                    valid = TRUE;
                }
            }
        }
    }

    pc->nextRow = pc->rowPos;
    pc->nextCol = pc->colPos;
    modifyDistance(map->q, &(map->trainerNodes[0]), map->time);
    map->terrain[pc->rowPos][pc->colPos].occupied = TRUE;
}

trainer_t createTrainer(map_t* map, int type)
{
    trainer_t trainer;
    int* movement = trainerMovement[type];

    int row = 0; 
    int col = 0;
    int direction = rand()%NUM_DIRECTIONS;

    if(type == PC)
    {
        trainer.numPokemon = 0; // The pokemon are cloned in the method movePC
    }
    else if(type == SENTRY)
    {
        while(map->terrain[row][col].type != DEADEND || map->terrain[row][col].occupied == 1)
        {
            row = rand()%map->rows;
            col = rand()%map->cols;
        }
    }
    else
    {
        while(movement[map->terrain[row][col].type] == INT_MAX || map->terrain[row][col].occupied == 1)
        {
            row = rand()%map->rows;
            col = rand()%map->cols;
        }
    }

    if(type == PACER)
    {
        if(movement[map->terrain[row-1][col].type] != INT_MAX)
        {
            direction = NORTH;
        }
        else if(movement[map->terrain[row+1][col].type] != INT_MAX)
        {
            direction = SOUTH;
        }
        else if(movement[map->terrain[row-1][col+1].type] != INT_MAX)
        {
            direction = EAST;
        }
        else if(movement[map->terrain[row][col-1].type] != INT_MAX)
        {
            direction = WEST;
        }
    }

    //pokemon_t* pokemon_team[TEAMSIZE];
    if(type != PC)
    {
        int n = 1;
        while(rand()%10 < 6 && n < 6) {n++;}
        trainer.numPokemon = n;

        for(int i = 0; i < n; i++)
        {
            trainer.pokemon_team[i] = createPokemon(map->distance);
        }
    }

    int* bag = (int*)malloc(NUM_ITEM_TYPES*sizeof(int));
    bag[POTIONS] = 3;
    bag[REVIVES] = 2;
    bag[POKEBALLS] = 5;

    trainer.bag = bag;
    //trainer.pokemon_team = pokemon_team;
    trainer.movement = movement;
    trainer.symbol = trainerSymbols[type];
    trainer.rowPos = row;
    trainer.colPos = col;
    trainer.type = type;
    map->terrain[row][col].occupied = 1;
    trainer.nextRow = row;
    trainer.nextCol = col;
    trainer.spawnTerrain = map->terrain[row][col];
    trainer.direction = direction;
    trainer.defeated = FALSE;
    trainer.color = COLOR_PAIR(COLOR_WHITE);

    return trainer;
} 

void spawnTrainers(map_t* map, int numTrainers)
{
    trainer_t* trainers = (trainer_t*)malloc((numTrainers+1)*sizeof(trainer_t));

    trainers[PC] = createTrainer(map, PC);
    for(int i = 1; i <= numTrainers; i++)
    {
        int type = (i-1)%(NUM_TRAINER_TYPES-1) + 1; //PC is a special case
        trainers[i] = createTrainer(map, type);
    }

    map->numTrainers = numTrainers;
    map->trainers = trainers;
}

void destroyTrainer(trainer_t* trainer)
{
    for(int i = 0; i < trainer->numPokemon; i++)
    {
        delete trainer->pokemon_team[i];
    }

    free(trainer->bag);
    free(trainer);
}

//Returns the completion time of the next move
int calcNextMove(map_t* map, trainer_t* trainer, int** distanceMap, int time)
{
    int curRow = trainer->rowPos;
    int curCol = trainer->colPos;

    int rowChange[NUM_DIRECTIONS] = {-1, 1, 0, 0}; //n, s, e, w
    int colChange[NUM_DIRECTIONS] = {0, 0, 1, -1};

    if(trainer->type == PACER)
    {
        int nextRow = curRow + rowChange[trainer->direction];
        int nextCol = curCol + colChange[trainer->direction];
        if(!(nextRow < 0 || nextCol < 0 || nextRow >= map->rows || nextCol >= map->cols 
        || map->terrain[nextRow][nextCol].occupied == TRUE || trainer->movement[map->terrain[nextRow][nextCol].type] == INT_MAX))
        {
            trainer->nextRow = nextRow;
            trainer->nextCol = nextCol;
        }
        else
        {
            if(trainer->direction == NORTH || trainer->direction == EAST) //N -> S, E -> W
            {
                trainer->direction++;
            }
            else                                                         //N <- S, E <- W
            {
                trainer->direction--;
            }
        }
    }
    else if(trainer->type == WANDERER)
    {
        int nextRow = curRow + rowChange[trainer->direction];
        int nextCol = curCol + colChange[trainer->direction];
        if(map->terrain[nextRow][nextCol].symbol == trainer->spawnTerrain.symbol && !(nextRow < 0 || nextCol < 0
        || nextRow >= map->rows || nextCol >= map->cols || map->terrain[nextRow][nextCol].occupied == TRUE 
        || trainer->movement[map->terrain[nextRow][nextCol].type] == INT_MAX))
        {
            trainer->nextRow = nextRow;
            trainer->nextCol = nextCol;
        }
        else
        {
            trainer->direction = rand()%NUM_DIRECTIONS;
        }
    }
    else if(trainer->type == SENTRY)
    {
        trainer->nextRow = curRow;
        trainer->nextCol = curCol;
    }
    else if(trainer->type == EXPLORER)
    {
        int nextRow = curRow + rowChange[trainer->direction];
        int nextCol = curCol + colChange[trainer->direction];
        if(!(nextRow < 0 || nextCol < 0 || nextRow >= map->rows || nextCol >= map->cols 
        || map->terrain[nextRow][nextCol].occupied == TRUE || trainer->movement[map->terrain[nextRow][nextCol].type] == INT_MAX))
        {
            trainer->nextRow = nextRow;
            trainer->nextCol = nextCol;
        }
        else
        {
            trainer->direction = rand()%NUM_DIRECTIONS;
        }
    }
    else //Rival and Hiker gradient descent method
    {
        int distance = distanceMap[curRow][curCol];
        for(int i = -1; i <= 1; i++)
        {
            for(int j = -1; j <= 1; j++)
            {
                if(!(trainer->defeated) && distanceMap[curRow+i][curCol+j] < distance)
                {
                    distance = distanceMap[curRow+i][curCol+j];
                    trainer->nextRow = curRow+i;
                    trainer->nextCol = curCol+j;
                }
                else if(trainer->defeated && distanceMap[curRow+i][curCol+j] > distance
                && distanceMap[curRow+i][curCol+j] - distanceMap[curRow][curCol] < 50)
                // Trainers will get trapped in buildings without the <50
                {
                    distance = distanceMap[curRow+i][curCol+j];
                    trainer->nextRow = curRow+i;
                    trainer->nextCol = curCol+j;
                }
                //Introduces some randomness to avoid getting stuck in local minima
                else if(distanceMap[curRow+i][curCol+j] == distance && rand()%2)
                {
                    trainer->nextRow = curRow+i;
                    trainer->nextCol = curCol+j;
                }
            }
        }
    }

    int finishTime = time + (trainer->movement[map->terrain[trainer->nextRow][trainer->nextCol].type]);
    return finishTime;
}

void printMap(map_t* map)
{
    erase();
    printw("\n");

    for(int i = 0; i < ROWS; i++)
    {
        for(int j = 0; j < COLS; j++)
        {
            if(map->terrain[i][j].occupied == 0)
            {
                attron(map->terrain[i][j].color);
                printw("%c", map->terrain[i][j].symbol);
                attroff(map->terrain[i][j].color);
            }
            else
            {
                for(int k = 0; k <= map->numTrainers; k++)
                {
                    if(map->trainers[k].rowPos == i && map->trainers[k].colPos == j)
                    {
                        attron(map->trainers[k].color);
                        printw("%c", map->trainers[k].symbol);
                        attroff(map->trainers[k].color);
                    }
                }
            }
        }
        printw("\n");
    }
    refresh();
}

void printTrainers(map_t* map, int startIndex)
{
    erase();
    printw("  Trainer Information - Press t to escape \n");
    printw("    Press UP ARROW to scroll upward   \n");

    for(int i = startIndex; i <= min(map->numTrainers, startIndex+ROWS); i++)
    {
        trainer_t pc = map->trainers[0];
        trainer_t trainer = map->trainers[i];

        int rowDiff = pc.rowPos - trainer.rowPos;
        char rowDir = 'N';
        if(rowDiff < 0)
        {
            rowDir = 'S';
        }

        int colDiff = pc.colPos - trainer.colPos;

        char colDir = 'E';
        if(colDiff > 0)
        {
            colDir = 'W';
        }

        printw("ID: %2d, Type: %c, Distance from PC: %2d%c %2d%c\n", 
        i, trainer.symbol, abs(rowDiff), rowDir, abs(colDiff), colDir);
    }

    printw("  Press DOWN ARROW to scroll downward  \n");
    refresh();
}

pokemon_t* swapPokemon(trainer_t* pc)
{
    int numFaintedPokemon = 0;
    move(9, 0);
    clrtobot();


    for(int i = 0; i < pc->numPokemon; i++)
    {
        pokemon_t* pokemon = pc->pokemon_team[i];
        if(pokemon->hp <= 0)
        {
            numFaintedPokemon++;
            attron(COLOR_PAIR(COLOR_GRAY));
            printw("LVL %d %s, HP: %d/%d \n", pokemon->level, pokemon->name, pokemon->hp, pokemon->stats[HP]);
            attroff(COLOR_PAIR(COLOR_GRAY));
        }
        else
        {
            printw("LVL %d %s, HP: %d/%d \n", pokemon->level, pokemon->name, pokemon->hp, pokemon->stats[HP]);
        }
    }

    if(numFaintedPokemon == pc->numPokemon)
    {
        pc->defeated = TRUE;
        pc->color = COLOR_PAIR(COLOR_GRAY);
        return nullptr;
    }

    int index = 0;
    int input = -1;
    while(input != 10 || pc->pokemon_team[index]->hp <= 0)
    {
        mvprintw(9+index, 40, "<----");
        refresh();
        input = getch();
        mvprintw(9+index, 40, "     ");

        if(input == KEY_DOWN && index < pc->numPokemon - 1) {index++;}
        else if(input == KEY_UP && index > 0) {index--;}
    }

    return pc->pokemon_team[index];
}

int revivePokemon(trainer_t* pc)
{
    move(9, 0);
    clrtobot();
    
    std::vector<pokemon_t*> faintedPals;
    for(int i = 0; i < pc->numPokemon; i++)
    {
        if(pc->pokemon_team[i]->hp <= 0)
        {
            faintedPals.push_back(pc->pokemon_team[i]);
            printw("Level %d %s\n", pc->pokemon_team[i]->level, pc->pokemon_team[i]->name);
        }
    }

    if(faintedPals.size() == 0)
    {
        return -1;
    }

    printw("BACK");

    int input = -1;
    int index = 0;
    while(input != 10)
    {
        mvprintw(9+index, 30, "<----");
        refresh();
        input = getch();
        mvprintw(9+index, 30, "     ");

        if(input == KEY_UP && index > 0) {index--;}
        else if(input == KEY_DOWN && index < (int)faintedPals.size()) {index++;}
    }

    if(index == (int)faintedPals.size())
    {
        return -1;
    }
    else
    {
        faintedPals[index]->hp = faintedPals[index]->stats[HP]/2;
        clrtobot();
        mvprintw(9, 0, "%s has been revived with %d HP.", faintedPals[index]->name, faintedPals[index]->hp);
        refresh();
        getch();
        return 1;
    }
}

int enterBag(trainer_t* pc)
{
    clrtobot();
    mvprintw(9, 0, "Potions: %d\n", pc->bag[POTIONS]);
    printw("Revives: %d\n", pc->bag[REVIVES]);
    printw("Pokeballs: %d\n", pc->bag[POKEBALLS]);
    printw("BACK");

    int input = - 1;
    int index = 0;
    while(input != 10)
    {
        mvprintw(9 + index, (int)strlen("Pokeballs: ")+3, "<----");
        refresh();
        input = getch();
        mvprintw(9 + index, (int)strlen("Pokeballs: ")+3, "     ");

        if(input == KEY_DOWN && index < NUM_ITEM_TYPES)
        {
            index++;
        }
        else if(input == KEY_UP && index > 0)
        {
            index--;
        }
    }
    
    if (index == NUM_ITEM_TYPES)
    {
        return -1;
    }
    else
    {
        if(pc->bag[index] <= 0)
        {
            mvprintw(14, 0, "None of those available.");
            refresh();
            getch();
            return -1;
        }
        pc->bag[index]--;
        return index;
    }
}

pokemon_t* getAlivePokemon(trainer_t* trainer)
{
    for(int i = 0; i < trainer->numPokemon; i++)
    {
        if(trainer->pokemon_team[i]->hp > 0)
        {
            return trainer->pokemon_team[i];
        }
    }

    trainer->defeated = TRUE;
    trainer->color = COLOR_PAIR(COLOR_GRAY);
    return nullptr;
}

void printBattleWindow(pokemon_t* trainer_pokemon, pokemon_t* pc_pokemon)
{
    mvprintw(3, 0,"LVL %d %s, HP: %d      LVL %d %s HP: %d \n", 
    pc_pokemon->level, pc_pokemon->name, pc_pokemon->hp, trainer_pokemon->level, trainer_pokemon->name, trainer_pokemon->hp);
    printw("------------------------------------\n");
    printw("               |                 \n");
    printw("      %c        |        %c       \n", pc_pokemon->name[0], trainer_pokemon->name[0]);
    printw("               |                 \n");
    printw("------------------------------------\n");
}

enum battlechoices {
    FIGHT,
    BAG,
    POKEMON,
    RUN
};

int promptBattleChoices()
{
    move(9, 0);
    clrtobot();
    printw("FIGHT           BAG \n");
    printw("POKEMON         RUN \n");

    int input = -1;
    int x = 0;
    int y = 0;
    while(input != 10)
    {
        mvprintw(9+y, 8+x*12, "<----");
        refresh();
        input = getch();
        mvprintw(9+y, 8+x*12, "     ");

        if(input == KEY_UP && y == 1)         {y--;}
        else if(input == KEY_RIGHT && x == 0) {x++;}
        else if(input == KEY_LEFT && x == 1)  {x--;}
        else if(input == KEY_DOWN && y == 0)  {y++;}
    }

    if(x == 0 && y == 0)
    {
        return FIGHT;
    }
    else if(x == 1 && y == 0)
    {
        return BAG;
    }
    else if(x == 0 && y == 1)
    {
        return POKEMON;
    }
    else
    {
        return RUN;
    }
}

int promptMoveChoice(pokemon_t* trainer_pokemon, pokemon_t* pc_pokemon)
{
    char* move0 = moves[pc_pokemon->moves[0]].identifier;
    char* move1 = moves[pc_pokemon->moves[1]].identifier;

    clrtobot();
    mvprintw(9, 0, "%s       %s       BACK", move0, move1);
    int cursorPos[3] = {int(strlen(move0)) + 1, int(strlen(move0)) + int(strlen(move1)) + 8, int(strlen(move0)) + int(strlen(move1)) + 19};

    int index = 0;
    int input = -1;
    while(input != 10)
    {
        mvprintw(9, cursorPos[index], "<----");
        refresh();
        input = getch();
        mvprintw(9, cursorPos[index], "     ");

        if(input == KEY_RIGHT && index < 2) {index++;}
        else if(input == KEY_LEFT && index > 0) {index--;}
    }

    return index;
}

void battle(trainer_t* pc, pokemon_t* pc_pokemon, pokemon_t* trainer_pokemon, bool isWild)
{
    while(trainer_pokemon->hp > 0 && pc_pokemon->hp > 0)  
    {
        printBattleWindow(trainer_pokemon, pc_pokemon);
        int action = promptBattleChoices();

        if(action == FIGHT)
        {
            int pc_move = promptMoveChoice(trainer_pokemon, pc_pokemon);
            int trainer_move = rand()%2;

            if(pc_move == 2)
            {
                continue;
            }
            else
            {
                if(moves[pc_pokemon->moves[pc_move]].priority > moves[trainer_pokemon->moves[trainer_move]].priority 
                || (moves[pc_pokemon->moves[pc_move]].priority == moves[trainer_pokemon->moves[trainer_move]].priority && pc_pokemon->stats[SPEED] > trainer_pokemon->stats[SPEED])
                || (moves[pc_pokemon->moves[pc_move]].priority == moves[trainer_pokemon->moves[trainer_move]].priority && pc_pokemon->stats[SPEED] > trainer_pokemon->stats[SPEED] && rand()%2 == 0))
                {
                    move(9,0);
                    clrtobot();
                    printw(attack(pc_pokemon, pc_move, trainer_pokemon).c_str());
                    refresh();
                    getch();

                    if(trainer_pokemon->hp < 0)
                    {
                        move(9,0);
                        clrtobot();
                        printw("%s has fainted.", trainer_pokemon->name);
                        refresh();
                        getch();
                        break;
                    }

                    move(9,0);
                    clrtobot();
                    printw(attack(trainer_pokemon, trainer_move, pc_pokemon).c_str());
                    refresh();
                    getch();

                    if(pc_pokemon->hp < 0)
                    {
                        move(9,0);
                        clrtobot();
                        printw("%s has fainted.", pc_pokemon->name);
                        refresh();
                        getch();
                    }
                }
                else
                {
                    move(9,0);
                    clrtobot();
                    printw(attack(trainer_pokemon, trainer_move, pc_pokemon).c_str());
                    refresh();
                    getch();

                    if(pc_pokemon->hp < 0)
                    {
                        move(9,0);
                        clrtobot();
                        printw("%s has fainted.", pc_pokemon->name);
                        refresh();
                        getch();
                    }
                    
                    move(9,0);
                    clrtobot();
                    printw(attack(pc_pokemon, pc_move, trainer_pokemon).c_str());
                    refresh();
                    getch();

                    if(pc_pokemon->hp < 0)
                    {
                        move(9,0);
                        clrtobot();
                        printw("%s has fainted.", trainer_pokemon->name);
                        refresh();
                        getch();
                        break;
                    }
                }
            }
        }
        else
        {
            if(action == RUN)
            {
                if(isWild)
                {
                    int oddsEscape = 100 - trainer_pokemon->level;
                    if(rand()%100 < oddsEscape)
                    {
                        printw("Got away safely.");
                        refresh();
                        getch();
                        break;
                    }
                    else
                    {
                        printw("Tried to run away but failed.");
                        refresh();
                        getch();
                    }
                }
                else
                {
                    mvprintw(9, 0, "\n\n");
                    mvprintw(9, 0, "Cannot run away from a trainer!");
                    refresh();
                    getch();
                    continue;
                }
            }
            else if(action == POKEMON)
            {
                pc_pokemon = swapPokemon(pc);
                printBattleWindow(trainer_pokemon, pc_pokemon);
                move(9, 0);
                clrtobot();
                printw("%s was sent out.", pc_pokemon->name);
                refresh();
                getch();
            }
            else if(action == BAG)
            {
                move(9,0);
                clrtobot();

                int item = enterBag(pc);
                move(9,0);
                clrtobot();

                if(item == POKEBALLS)
                {
                    if(isWild == FALSE)
                    {
                        printw("Cannot capture another trainer's pokemon!");
                        refresh();
                        getch();
                        pc->bag[item]++;
                        continue;
                    }
                    else
                    {
                        if(pc->numPokemon < 6 && (rand()%3 == 0))
                        {
                            pc->pokemon_team[pc->numPokemon] = trainer_pokemon;
                            pc->numPokemon++;

                            printw("The wild %s has been captured.", trainer_pokemon->name);
                            refresh();
                            getch();
                            break;
                        }
                        else
                        {
                            printw("The wild %s has has escaped your pokeball.", trainer_pokemon->name);
                            refresh();
                            getch();
                        }
                    }
                }
                else if(item == POTIONS)
                {
                    int healthGained = std::min(20, pc_pokemon->stats[HP] - pc_pokemon->hp);
                    printw("Used a potion on %s.\n", pc_pokemon->name);
                    printw("%s gained %d health.", pc_pokemon->name, healthGained);
                    pc_pokemon->hp += healthGained;
                    refresh();
                    getch();
                }
                else if(item == REVIVES)
                {
                    int pokemonIndex = revivePokemon(pc);
                    if(pokemonIndex == -1)
                    {
                        printw("No pokemon to revive.");
                        refresh();
                        getch();
                        pc->bag[item]++;
                        continue;
                    }
                }
                else {continue;} //BACK was selected
            }

            int trainer_move = rand()%2;
            move(9,0);
            clrtobot();
            printw(attack(trainer_pokemon, trainer_move, pc_pokemon).c_str());
            refresh();
            getch();

            if(pc_pokemon->hp < 0)
            {
                move(9,0);
                clrtobot();
                printw("%s has fainted.", pc_pokemon->name);
                refresh();
                getch();
            }
        }
    }
}

void battlePokemon(trainer_t* pc, pokemon_t* pokemon)
{
    bool isWild = true;
    pokemon_t* pc_pokemon = getAlivePokemon(pc);

    erase();
    printw("A wild %s appears!", pokemon->name);

    while(pc->defeated == false)
    {
        battle(pc, pc_pokemon, pokemon, isWild);
        
        if(pc_pokemon->hp <= 0)
        {
            pc_pokemon = swapPokemon(pc);
        }
        else if(pokemon->hp <= 0)
        {
            delete pokemon;
            break;
        }
        else
        {
            break;
        }
    }

    return;
}

int healPokemon(trainer_t* pc)
{
    erase();
    for(int i = 0; i < pc->numPokemon; i++)
    {
        pokemon_t* pokemon = pc->pokemon_team[i];
        printw("LVL %d %s, HP: %d/%d \n", pokemon->level, pokemon->name, pokemon->hp, pokemon->stats[HP]);
    }
    printw("BACK");

    int input = -1;
    int index = 0;
    while(input != 10)
    {
        mvprintw(index, 40, "<----");
        refresh();
        input = getch();
        mvprintw(index, 40, "     ");

        if(input == KEY_UP && index > 0) {index--;}
        else if(input == KEY_DOWN && index < pc->numPokemon) {index++;}
    }

    if(index == pc->numPokemon || (pc->pokemon_team[index]->hp == pc->pokemon_team[index]->stats[HP]))
    {
        return FALSE;
    }
    else
    {
        int healthGained = std::min(20, pc->pokemon_team[index]->stats[HP] - pc->pokemon_team[index]->hp);
        pc->pokemon_team[index]->hp += healthGained;

        erase();
        printw("Used a potion on %s.\n", pc->pokemon_team[index]->name);
        printw("%s gained %d health.", pc->pokemon_team[index]->name, healthGained);
        getch();
        return TRUE;
    }
}

void battleTrainer(trainer_t* pc, trainer_t* trainer)
{
    erase();

    bool isWild = false;
    pokemon_t* pc_pokemon = getAlivePokemon(pc);

    for(int i = 0; i < trainer->numPokemon; i++)
    { 
        pokemon_t* trainer_pokemon = trainer->pokemon_team[i];

        erase();
        printw("Initiated battle with %c! \n", trainer->symbol);
        printw("%c sent out %s! {%d/%d}\n\n", trainer->symbol, trainer_pokemon->name, trainer->numPokemon-i, trainer->numPokemon);
        battle(pc, pc_pokemon, trainer_pokemon, isWild);

        if(pc_pokemon->hp <= 0)
        {
            pc_pokemon = swapPokemon(pc);
            if(pc->defeated == TRUE)
            {
                return;
            }
            else
            {
                i--;
            }
        }
    }

    trainer->defeated = TRUE;
    trainer->color = COLOR_PAIR(COLOR_GRAY);
    return;
}

#define BATTLE 2
#define GATE 3
int validNextMove(map_t* map, trainer_t* trainer)
{
    int nextRow = trainer->nextRow;
    int nextCol = trainer->nextCol;

    if(nextRow >= 0 && nextCol >= 0 && nextRow < map->rows && nextCol < map->cols 
    && (map->terrain[nextRow][nextCol].occupied == FALSE || (nextRow == trainer->rowPos && nextCol == trainer->colPos))
    && trainer->movement[map->terrain[nextRow][nextCol].type] != INT_MAX)
    {
        if(nextRow == 0 || nextCol == 0 || nextRow == map->rows-1 || nextCol == map->cols-1)
        {
            return GATE;
        }

        return TRUE;
    }
    else if(map->terrain[nextRow][nextCol].occupied)
    {
        trainer_t* pc = &(map->trainers[PC]);
        trainer_t* defender;
        for(int i = 0; i <= map->numTrainers; i++)
        {
            if(map->trainers[i].colPos == nextCol && map->trainers[i].rowPos == nextRow)
            {
                defender = &(map->trainers[i]);
            }
        }

        if(trainer->type == PC && (defender->defeated == FALSE))
        {
            battleTrainer(pc, defender);
            return BATTLE;
        }
        else if(defender->type == PC && (trainer->defeated == FALSE))
        {
            battleTrainer(pc, trainer);
            return BATTLE;
        }
        else
        {
            return FALSE;
        }
    }
    else
    {
        return FALSE;
    }
}

void moveTrainer(map_t* map, trainer_t* trainer)
{
    int validMove = validNextMove(map, trainer);
    if(validMove == FALSE || validMove == BATTLE)
    {
        return;
    }
    else if(validMove == GATE)
    {
        map->terrain[trainer->rowPos][trainer->colPos].occupied = FALSE;
    }
    else
    {
        map->terrain[trainer->rowPos][trainer->colPos].occupied = FALSE;
        trainer->rowPos = trainer->nextRow;
        trainer->colPos = trainer->nextCol;
        map->terrain[trainer->rowPos][trainer->colPos].occupied = TRUE;
    }

    return;
}

// Hi