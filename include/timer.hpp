//
// Created by Giovanni Bollati on 09/11/23.
//

#ifndef GL_TEST_PROJECT_TIMER_HPP
#define GL_TEST_PROJECT_TIMER_HPP

#include <iostream>
#include <SDL.h>

class Timer {
private:
    int id;
    Uint64 ms;
public:
    Timer() : ms(SDL_GetTicks64()) {
        static int _id = 0;
        id = _id++;
    }
    int time() { return SDL_GetTicks64() - ms; }

    void log() {
        std::cout << std::endl << "timer " << id;
        if(ms > 0) std::cout << " | time elapsed: " << time();
        std::cout << std::endl;
    }
};

#endif //GL_TEST_PROJECT_TIMER_HPP
