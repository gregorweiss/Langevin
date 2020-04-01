//
// Created by gregor on 14.07.19.
//

#include "ArgParse.h"

ArgParse::ArgParse(int argc, char *argv[]) :
        _argc(argc),
        _argv(new char *[argc * sizeof(char)]),
        iMap() {

    for (int i = 0; i < argc; i++) {
        _argv[i] = new char[strlen(argv[i]) + 1];
        strncpy(_argv[i], argv[i], strlen(argv[i]));
        _argv[i][strlen(argv[i])] = '\0';
    }
}
