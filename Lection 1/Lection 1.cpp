// Lection 1.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <algorithm>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <windows.h>
#include <intrin.h>
using namespace std;

//NTSTATUS NtQueryTimerResolution(OUT PULONG MinimumResolution, OUT PULONG MaximumResolution, OUT PULONG CurrentResolution);
//NTSTATUS  TimesResolution(PULONG pMinimumResolution, OUT PULONG pMaximumResolution, OUT PULONG pCurrentResolution)
//{
//    NTSTATUS status = NtQueryTimerResolution(pMinimumResolution, pMaximumResolution, pCurrentResolution);
//    return status;
//}
#if 0
typedef INT(WINAPI* tNtQueryTimerResolution)(PULONG, PULONG, PULONG);
NTSTATUS  TimesResolution(PULONG pMinimumResolution, OUT PULONG pMaximumResolution, OUT PULONG pCurrentResolution)
{    
    NTSTATUS status;
    HMODULE hDLL = GetModuleHandleW(L"ntdll.dll");
    if (hDLL) {
        tNtQueryTimerResolution NtQueryTimerResolution =
            (tNtQueryTimerResolution)GetProcAddress(hDLL, "NtQueryTimerResolution");
        if (NtQueryTimerResolution) {
            status = NtQueryTimerResolution(pMinimumResolution, pMaximumResolution, pCurrentResolution);
            
        }
    }
    return status;
}

double GetAccuracyClock() {
    clock_t Start, Finish;
    Start = clock();
    do { Finish = clock(); } while (Finish - Start == 0);
    double Accurancy = (double)(Finish - Start) / CLOCKS_PER_SEC;
    return Accurancy;
}	

#include <chrono>
using namespace std::chrono;
double GetChronoAccuracy() {
    auto start = high_resolution_clock::now();
    auto finish = high_resolution_clock::now();
    do {
        finish = high_resolution_clock::now();
    } while ((finish - start).count() == 0);
    auto dif = std::chrono::duration_cast<std::chrono::nanoseconds>
        (finish - start).count();
    return double(dif)/ 1e9;				//100 ns = 0.0001 ms
}

#include <omp.h>
double GetOMPAccuracy() {
    double Start, Finish;
    Start = omp_get_wtime();
    do {
        Finish = omp_get_wtime();
    } while (Finish - Start == 0);
    return (Finish - Start);
}





#pragma region task1
#define N   (1 << 16)
void bubbleSort(int arr[], int n) {
    for (int i = 0; i < n - 1; i++) {
        for (int j = 0; j < n - i - 1; j++) {
            if (arr[j] > arr[j + 1]) {
                int temp = arr[j]; arr[j] = arr[j + 1]; arr[j + 1] = temp;
            }
        }
    }
}
#pragma endregion task1
#endif
#pragma region task2
#define M   2048
void multiplyMatrices(float A[][M], float B[][M], 
    float result[][M], int size) {
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            result[i][j] = 0.0;
            for (int k = 0; k < size; k++) {
                result[i][j] += A[i][k] * B[k][j];
            }
        }
    }	
}

void multiplyMatrices1(float A[][M], float B[][M], float result[][M], int size) {
    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++)  result[i][j] = 0.0;
    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++)
            for (int k = 0; k < size; k++)
                result[i][j] += A[i][k] * B[k][j];
}

void multiplyMatrices2(float A[][M], float B[][M], float result[][M], int size) {
    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++)  result[i][j] = 0.0;
    for (int i = 0; i < size; i++)
        for (int k = 0; k < size; k++)
            for (int j = 0; j < size; j++)
                result[i][j] += A[i][k] * B[k][j];
}
#pragma endregion task2
#pragma region task3
#include <cmath>
void std_round(float* dest, float* src, size_t n) {
    size_t i;
    for (i = 0; i < n; i++)	dest[i] = std::round(src[i]);
}

void my_round(float* dest, float* src, size_t n) {
    size_t i;
    for (i = 0; i < n; i++) {
        if (src[i] < 0)
            dest[i] = (float)(int64_t)(src[i] - 0.5);
        else
            dest[i] = (float)(int64_t)(src[i] + 0.5);
    }
}

#include <intrin.h>
void avx_round(float* dest, float* src, size_t n) {
    __m256* dest256 = (__m256*)dest;
    __m256* src256 = (__m256*)src;
    size_t i;
    for (i = 0; i < n / 8; i++)
        dest256[i] = _mm256_round_ps(src256[i], _MM_FROUND_TO_NEAREST_INT);
}

#include <omp.h>
void my_roundOMP(float* dest, float* src, size_t n) {
#pragma omp parallel for
    for (int i = 0; i < n; i++) {
        if (src[i] < 0)
            dest[i] = (float)(int64_t)(src[i] - 0.5);
        else
            dest[i] = (float)(int64_t)(src[i] + 0.5);
    }
}


#pragma endregion task3



//#define N   (1024)
//int etalon[N];
//int esort[N], unsort[N], cur[N];
float m1[M][M], m2[M][M], m3[M][M], m4 [M][M];
int main()
{
    clock_t t1, t2, t3;
#if 0
    for (int i = 0; i < N; ++i)
        etalon[i] = rand();
    std::copy(etalon, etalon + N, cur);
    std::copy(etalon, etalon + N, esort);
    std::sort(esort, esort + N);
    for (int i = N - 1, j = 0; i >= 0; i--)
    {
        unsort[j++] = esort[i];
    }
    
    
    

    t1 = clock();
    bubbleSort(esort, N);
    t1 = clock() - t1;

    t3 = clock();
    bubbleSort(cur, N);
    t3 = clock() - t3;

    t2 = clock();
    bubbleSort(unsort, N);
    t2 = clock() - t2;


    
    printf("n = %d\n", N);
    printf("Waiting t1 < t3 < t2\n");
    printf("t1 = %d t3 = %d t2 = %d \n",
        t1, t3, t2);
    
    printf("esort [N - 1] = %d unsort [N - 1] = %d cur [N - 1] = %d\n",
        esort[N - 1], unsort[N - 1], cur[N - 1]);
#endif

#if 1
    for (int i = 0; i < M; ++i)
    {
        for (int j = 0; j < M; j++)
        {
            m1[i][j] = rand() % 10 + 0.f;
            m2[i][j] = rand() % 10 + 0.f;

        }
    }
    printf("M = %d\n", M);
#endif
#if 1
    t1 = clock();
    multiplyMatrices(m1, m2, m3, M);
    t1 = clock() - t1;
    printf("m3[0][0] = %g m3[M-1][M-1] = %g\n", m3[0][0], m3[M - 1][M - 1]);
    t2 = clock();
    multiplyMatrices1(m1, m2, m4, M);
    t2 = clock() - t2;
    printf("m4[0][0] = %g m4[M-1][M-1] = %g\n", m4[0][0], m4[M - 1][M - 1]);

    t3 = clock();
    multiplyMatrices2(m1, m2, m4, M);
    t3 = clock() - t3;
    printf("m4[0][0] = %g m4[M-1][M-1] = %g\n", m4[0][0], m4[M - 1][M - 1]);
    printf("t1 = %d t2 = %d t3 = %d\n", t1, t2, t3);
#endif
#if 0
    float* pm1 = (float*)m1, *pm2 = (float*)m2, * pm3 = (float*)m3;
    for (int i = 0; i < M * M; ++i)
    {
        pm1[i] = rand()  / (rand() + 1.f);
    }
    bool b = true;
    t1 = clock();
    std_round(pm2, pm1, M * M);
    t1 = clock() - t1;
    t2 = clock();
    my_roundOMP(pm3, pm1, M * M);
    // my_round(pm3, pm1, M * M);
    t2 = clock() - t2;
    for (int i = 0; i < M * M; ++i) {
        if (pm2[i] != pm3[i])
        {
            b = false;
            printf("std_round and my_round: i = %d pm1 [i] = %g pm2 [i] = %g pm3 [i] = %g\n", 
                i, pm1[i], pm2[i], pm3[i]);
            break;
        }
    }
    printf("Compare std_round and my_round ? %s\n", b ? "OK" : "ERROR");
    b = true;
    t3 = clock();
    avx_round(pm3, pm1, M * M);
    t3 = clock() - t3;
    for (int i = 0; i < M * M; ++i) {
        if (pm2[i] != pm3[i])
        {
            b = false;
            printf("std_round and avx_round: i = %d pm1 [i] = %g pm2 [i] = %g pm3 [i] = %g\n",
                i, pm1[i], pm2[i], pm3[i]);
            break;
        }
    }
    printf("Compare std_round and avx_round ? %s\n", b ? "OK" : "ERROR");
    printf("t1 = %d t2 = %d t3 = %d\n", t1, t2, t3);
#endif

    //b = true;
    //t3 = clock();
    //my_roundOMP(pm3, pm1, M * M);
    ////my_round(pm3, pm1, M * M);
    //t3 = clock() - t3;
    //for (int i = 0; i < M * M; ++i) {
    //    if (pm2[i] != pm3[i])
    //    {
    //        b = false;
    //        printf("std_round and my_roundOMP: i = %d pm1 [i] = %g pm2 [i] = %g pm3 [i] = %g\n",
    //            i, pm1[i], pm2[i], pm3[i]);
    //        break;
    //    }
    //}
    //printf("Compare std_round and my_roundOMP ? %s\n", b ? "OK" : "ERROR");
    //printf("t3 = %d\n", t3);
#if 0 
    ULONG MinimumResolution, MaximumResolution, CurrentResolution;
    NTSTATUS  status = TimesResolution(&MinimumResolution, &MaximumResolution, &CurrentResolution);
    printf("MinimumResolution = %I64d MaximumResolution = %I64d CurrentResolution = %I64d\n",
        MinimumResolution, MaximumResolution, CurrentResolution);

    double AccuracyClock = GetAccuracyClock();
    printf("AccuracyClock = %lg\n", AccuracyClock);         // AccuracyClock = 0.001
    double AccuracyChrono = GetChronoAccuracy();
    printf("AccuracyChrono = %lg\n", AccuracyChrono);       // AccuracyChrono = 1e-07
    double AccuracyOMP = GetOMPAccuracy();
    printf("AccuracyOMP = %lg\n", AccuracyOMP);             // AccuracyOMP = 9.98843e-08
#endif
            
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
