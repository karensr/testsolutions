// ConsoleApplication1.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include "stdlib.h"

/*
#define MAX_N 100000
#define MAX_VAL 1000000

#define sullb (sizeof(unsigned long long) * 8)

#define numOfULL ((MAX_VAL + sullb - 1) / sullb)

unsigned long long dataAbit[numOfULL];

int solution(int A[], int N) {
    // write your code in C99 (gcc 6.2.0)

    int i, j;
    unsigned long long tmp = 0;

    dataAbit[0] = 1;

    for (i = 0; i < N; i++) {
        if (A[i] > 0) {
            int element = A[i] / sullb;
            int bit = A[i] % sullb;
            dataAbit[element] |= ((unsigned long long)1) << bit;
        }
    }
    
    for (i = 0; i < numOfULL; i++) {
        if (dataAbit[i] != 0xFFFFFFFFFFFFFFFF) {
            tmp = dataAbit[i];
            break;
        }
    }

    for (j = 0; j < sullb; j++) {
        if ((tmp & (unsigned long long)1) == 0) {
            return (i * sullb + j);
        }
        tmp = tmp >> 1;
    }
    
    return MAX_VAL;
}
*/
/*
#define sintb (sizeof(int) * 8)

int solution(int N) {
    int ff = -1;
    int gap = 0;
    int tmpgap = 0;
    if (N < 5)
        return 0;

    for (int i = 0; i < sintb; i++)
    {
        if (N & 1) {
            if (ff == -1) {
                ff = i;
            }
            else if ((i - ff) > 1) {
                tmpgap = i - ff - 1;
                if (tmpgap > gap) {
                    gap = tmpgap;
                }
            }
            ff = i;
        }
        N = N >> 1;
    }
    return gap;
}
*/
/*
void myMemCpy(int* dest, int* src, int size)
{
    int* csrc = (int*)src;
    int* cdest = (int*)dest;

    for (int i = 0; i < size; i++)
        cdest[i] = csrc[i];
}

#define MAX_N 100

struct Results {
    int* A;
    int N; // Length of the array
};

int tmpArray[MAX_N];

struct Results solution(int A[], int N, int K) {
    static struct Results result;

    if (A == NULL || N == 0) {
        result.A = NULL;
        result.N = 0;
        return result;
    }

    if (N > MAX_N) {
        N = MAX_N;
    }
    K = K % N;
    myMemCpy(&tmpArray[0], &A[N - K], K);
    myMemCpy(&tmpArray[K], &A[0], N - K);
    // write your code in C99 (gcc 6.2.0)

    result.A = tmpArray;
    result.N = N;
    return result;
}
*/
/*
#define MAX_N 1000000
#define MAX_VAL 1000000000

#define sullb (sizeof(unsigned long long) * 8)

#define numOfULL ((MAX_VAL + sullb) / sullb)

unsigned long long dataAbit[numOfULL];

int solution(int A[], int N) {
    // write your code in C99 (gcc 6.2.0)

    int i, j;
    unsigned long long tmp = 0;

    for (i = 0; i < N; i++) {
        if (A[i] <= MAX_VAL && A[i] > 1) {
            int element = A[i] / sullb;
            int bit = A[i] % sullb;
            dataAbit[element] ^= (((unsigned long long)1) << bit);
        }
    }

    for (i = 0; i < numOfULL; i++) {
        if (dataAbit[i] != 0) {
            tmp = dataAbit[i];
            break;
        }
    }

    for (j = 0; j < sullb; j++) {
        if ((tmp & (unsigned long long)1) == 1) {
            return (i * sullb + j);
        }
        tmp = tmp >> 1;
    }

    return 0;
}
*/
/*
int solution(int X, int Y, int D) {
    if (X >= Y) {
        return 0;
    }

    return ((Y - X + D - 1) / D);
    // write your code in C99 (gcc 6.2.0)
}
*/
/*
#define MAX_N 100000

#define sullb (sizeof(unsigned long long) * 8)

#define numOfULL ((MAX_N + sullb) / sullb)

unsigned long long dataAbit[numOfULL];

int solution(int A[], int N) {
    // write your code in C99 (gcc 6.2.0)

    int i, j;
    unsigned long long tmp = 0;
    
    dataAbit[0] = 1;

    for (i = 0; i < N; i++) {
        int element = A[i] / sullb;
        int bit = A[i] % sullb;
        dataAbit[element] |= (((unsigned long long)1) << bit);
    }

    for (i = 0; i < numOfULL; i++) {
        if (dataAbit[i] != 0xFFFFFFFFFFFFFFFF) {
            tmp = dataAbit[i];
            break;
        }
    }

    for (j = 0; j < sullb; j++) {
        if ((tmp & (unsigned long long)1) == 0) {
            return (i * sullb + j);
        }
        tmp = tmp >> 1;
    }

    return (MAX_N + 1);
}
*/
/*
int diffMod(int a, int b) {
    a -= b;
    return ((a >= 0) ? a : -a);
}

int solution(int A[], int N) {
    int fsum = 0;
    int ssum = 0;
    int diff = 0;
    int tmpdiff = 100000*1000; //keep maximum diff value

    for (int i = 0; i < N; i++) {
        ssum += A[i];
    }

    for (int i = 0; i < (N-1); i++) {
        ssum = ssum - A[i];
        fsum += A[i];
        diff = diffMod(fsum, ssum);
        if (diff == 0) {
            return 0;
        }
        
        if (diff < tmpdiff) {
            tmpdiff = diff;
        }
    }
    return tmpdiff;
}
*/
/*
#define MAX_N 100000

#define sullb (sizeof(unsigned long long) * 8)

#define numOfULL ((MAX_N + sullb) / sullb)

unsigned long long dataAbit[numOfULL];

int solution(int X, int A[], int N) {
    int i;
    int active = 0;
    int found = 0;
    dataAbit[0] = 1;
    for (i = 0; i < N; i++) {
        int element = A[i] / sullb;
        int bit = A[i] % sullb;
        dataAbit[element] |= (((unsigned long long)1) << bit);
        if (A[i] == X) {
            active |= 1;
        }
        if ((X - 1) == i) {
            active |= 2;
        }
        if (active == 3) {
            int xcheck = X;
            int j = 0;
            found = 1;
            while (xcheck > 63) {
                xcheck -= 64;
                if (dataAbit[j] != 0xFFFFFFFFFFFFFFFF) {
                    found = 0;
                    j++;
                    break;
                }
                j++;
            }
            if (found == 1) {
                unsigned long long cmpdata = (((unsigned long long)1) << (X + 1)) - 1;
                if (dataAbit[j] != cmpdata) {
                    found = 0;
                }
            }
        }
        if (found) {
            break;
        }
    }
    if (i == N) {
        i = -1;
    }
    return i;
}
*/
/*
struct Results {
    int* C;
    int L; // Length of the array
};

#define MAX_N 100000

int arrayN[MAX_N + 1];

struct Results solution(int N, int A[], int M) {
    struct Results result;
    // write your code in C99 (gcc 6.2.0)

    int maxNum = 0;
    int maxN = N + 1;
    int* pA = A;
    int* pD;
    int maxDone = 0;

    result.C = &arrayN[1];
    result.L = N;

    while (M--) {
        if (*pA == maxN) {
            maxDone = maxNum;
        }
        else {
            pD = &arrayN[*pA];
            if (*pD < maxDone) {
                *pD = maxDone + 1;
            }
            else {
                (*pD)++;
            }
            if (*pD > maxNum) {
                maxNum = *pD;
            }
        }
        pA++;
    }

    pD = arrayN;
    while (N--) {
        if (*(++pD) < maxDone) {
            *pD = maxDone;
        }
    }

    return result;
}
*/
/*
#define MAX_N 100000

#define sullb (sizeof(unsigned long long) * 8)

#define numOfULL ((MAX_N + sullb) / sullb)

unsigned long long dataAbit[numOfULL];

int solution(int A[], int N) {
    // write your code in C99 (gcc 6.2.0)

    int i;
    int j = 0;
    int found = 1;

    dataAbit[0] = 1;

    for (i = 0; i < N; i++) {
        if (A[i] > MAX_N) {
            continue;
        }
        int element = A[i] / sullb;
        int bit = A[i] % sullb;
        dataAbit[element] |= (((unsigned long long)1) << bit);
    }

    while (N > 63) {
        N -= 64;
        if (dataAbit[j] != 0xFFFFFFFFFFFFFFFF) {
            found = 0;
            j++;
            break;
        }
        j++;
    }

    if (found == 1) {
        unsigned long long cmpdata = (((unsigned long long)1) << (N + 1)) - 1;
        if (dataAbit[j] != cmpdata) {
            found = 0;
        }
    }

    if (found) {
        return 1;
    }

    return 0;
}
*/
/*
int solution(int A, int B, int K) {
    // write your code in C99 (gcc 6.2.0)
    int xa = A / K;
    int sum = 0;
    if (A >= K && (A % K) == 0) {
        xa--;
    }
    if (A == 0) {
        sum++;
    }
    int xb = B / K;
    sum += xb - xa;

    return sum;
}
*/
/*
struct Results {
    int* A;
    int M; // Length of the array
};


#define MAX_M 50000
#define MAX_N 100000

int arrayN[MAX_N];
int arrayM[MAX_M];
int map[16];
int mapN[20];

struct Results solution(char* S, int P[], int Q[], int M) {
    struct Results result;
    // write your code in C99 (gcc 6.2.0)
    int i;
    int j;
    int tmp;
    
    mapN['A' - 'A'] = 1;
    mapN['C' - 'A'] = 2;
    mapN['G' - 'A'] = 4;
    mapN['T' - 'A'] = 8;

    for (i = 0; i < 16; i++) {
        if (i & 1) {
            map[i] = 1;
        } else
        if (i & 2) {
            map[i] = 2;
        } else
        if (i & 4) {
            map[i] = 3;
        }
        else {
            map[i] = 4;
        }
    }

    i = 0;
    while(S[i]) {
        arrayN[i] = mapN[S[i] - 'A'];
        i++;
    }

    for (i = 0; i < M; i++) {
        tmp = 0;
        for (j = P[i]; j <= Q[i]; j++) {
            tmp |= arrayN[j];
            if (tmp & 1) {
                break;
            }
        }
        arrayM[i] = map[tmp];
    }

    result.A = arrayM;
    result.M = M;
    return result;
}
*/
/*
double getLastMin(int A[], int N) {
    double average = 10000;
    double tmpaverage = 10000;
    int sum = A[0];
    for (int i = 1; i < N; i++) {
        sum += A[i];
        tmpaverage = ((double)sum) / (i + 1);
        if (i == 1) {
            average = tmpaverage;
        }
        else {
            if (tmpaverage < average) {
                average = tmpaverage;
            }
            else {
                break;
            }
        }
    }
    return average;
}

int solution(int A[], int N) {
    double min = 10001;
    double lastmin;
    int minnum = 0;
    for (int i = 0; i < (N - 1); i++) {
        lastmin = getLastMin(&A[i],N-i);
        if (lastmin < min) {
            min = lastmin;
            minnum = i;
        }
    }
    return minnum;
}
*/
/*
int solution(int A[], int N) {
    unsigned int num0 = 0;
    unsigned int num1 = 0;
    unsigned int numchanged = 0;
    unsigned int i;
    unsigned int count = 0;
    
    if (N == 1) {
        return 0;
    }
    for (i = 0; i < (N - 1); i++) {
        if (A[i]) {
            num1++;
        }
        if (numchanged == 3) {
            continue;
        }
        if (A[i] != A[i + 1]) {
            numchanged++;
        }
    }
    if(A[N-1]) {
        num1++;
    }
    num0 = N - num1;

    if (num0 == 0 || num1 == 0) {
        return 0;
    }

    if (numchanged == 1) {
        if (A[0]) {
            return 0;
        }
        count = num0 * num1;
        if (count > 1000000000) {
            return -1;
        }
        return count;
    }

    for (i = 0; i < (N - 1); i++) {
        if (A[i]) {
            num1--;
        }
        else {
            count += num1;
            if (count > 1000000000) {
                return -1;
            }
        }
    }
    return count;
}
*/
/*
unsigned char map[2*1000000+1];

int solution(int A[], int N) {
    // write your code in C99 (gcc 6.2.0)
    int count = 0;
    for(int i = 0; i < N; i++) {
        if(map[A[i] + 1000000]) {
            continue;
        }
        else {
            map[A[i] + 1000000] = 1;
            count++;
        }
    }
    return count;
}
*/
/*
int solution(int A[], int N) {
    int maxn1 = -1000;
    int maxn2 = -1000;
    int maxn3 = -1000;
    int minn1 = 1000;
    int minn2 = 1000;
    int val1, val2;

    for (int i = 0; i < N; i++) {
        if (A[i] >= maxn1) {
            maxn3 = maxn2;
            maxn2 = maxn1;
            maxn1 = A[i];
        }
        else
            if (A[i] >= maxn2) {
                maxn3 = maxn2;
                maxn2 = A[i];
            }
            else
                if (A[i] >= maxn3) {
                    maxn3 = A[i];
                }

        if (A[i] <= minn1) {
            minn2 = minn1;
            minn1 = A[i];
        }
        else
            if (A[i] <= minn2) {
                minn2 = A[i];
            }
    }

    val1 = minn1 * minn2 * maxn1;
    val2 = maxn1 * maxn2 * maxn3;

    return (val1 > val2 ? val1 : val2);
}
*/
/*
int solution(int A[], int N) {
    int counter = 0;
    int i, j;
    for (i = 0; i < (N - 1); i++) {
        for (j = 1; j < (N - i); j++) {
            if (A[i] >= j - A[j + i]) {
                counter++;
            }
        }
        if (counter > 10000000) {
            return -1;
        }
    }
    return counter;
}
*/
/*
/////////////////////////////////////////////////////////
void swap(int* a, int* b) {
    int temp = *a;
    *a = *b;
    *b = temp;
}

void quicksort(int arr[], int l, int r) {
    if (l >= r) {
        return;
    }
    int pivot = arr[r];
    int cnt = l;
    for (int i = l; i <= r; i++) {
        if (arr[i] <= pivot) {
            swap(&arr[cnt], &arr[i]);
            cnt++;
        }
    }

    quicksort(arr, l, cnt - 2);
    quicksort(arr, cnt, r);
}


int solution(int A[], int N) {
    quicksort(A, 0, N - 1);
    for (int i = 0; i < (N - 2); i++) {
        if ((long long)A[i] + A[i + 1] > A[i + 2]) {
            return 1;
        }
    }
    return 0;
}
/////////////////////////////// fast
int compare(const void* a, const void* b)
{
    long long la = *(int*)a;
    long long lb = *(int*)b;
    if (la > lb)
        return 1;
    if (la < lb)
        return -1;
    return 0;

}

int solution(int A[], int N) {
    qsort(A, N, sizeof(int), compare);
    for (int i = 0; i < (N - 2); i++) {
        if ((long long)A[i] + A[i + 1] > A[i + 2]) {
            return 1;
        }
    }
    return 0;
}
*/
/*
/*
struct sNode
{
    char data;
    struct sNode* next;
};

void push(struct sNode** top_ref, int new_data) {
    struct sNode* new_node = (struct sNode*) malloc(sizeof(struct sNode));
    if (new_node == NULL)
        return;
    new_node->data = new_data;
    new_node->next = (*top_ref);
    (*top_ref) = new_node;
}

int pop(struct sNode** top_ref) {
    char res;
    struct sNode* top;
    if (*top_ref == NULL) {
        printf("Stack overflow n");
        getchar();
        exit(0);
    }
    else {
        top = *top_ref;
        res = top->data;
        *top_ref = top->next;
        free(top);
        return res;
    }
}

int isMatchingPair(char character1, char character2) {
    if (character1 == '(' && character2 == ')')
        return 1;
    else if (character1 == '{' && character2 == '}')
        return 1;
    else if (character1 == '[' && character2 == ']')
        return 1;
    else
        return 0;
}

int solution(char* S) {
    int i = 0;
    struct sNode* stack = NULL;
    while (S[i]) {
        if (S[i] == '{' || S[i] == '(' || S[i] == '[')
            push(&stack, S[i]);

        if (S[i] == '}' || S[i] == ')' || S[i] == ']') {
            if (stack == NULL)
                return 0;
            else if (!isMatchingPair(pop(&stack), S[i]))
                return 0;
        }
        i++;
    }

    if (stack == NULL)
        return 1;

    return 0;
}
*/
/*
int solution(int A[], int B[], int N) {
    int i;
    int count = N - 1;
    int eatDetected = 0;
    do {
        count -= eatDetected;
        eatDetected = 0;
        i = 0;
        while (i != count) {
            if (B[i] == 0) {
                i++;
            }
            else {
                break;
            }
        }
        for (; i < count; i++) {
            if (B[i] != B[i + 1]) {
                if (A[i] > A[i + 1]) {
                    A[i - eatDetected] = A[i];
                    B[i - eatDetected] = B[i];
                    i++;
                }
                else {
                    A[i - eatDetected] = A[i + 1];
                    B[i - eatDetected] = B[i + 1];
                    i++;
                }
                eatDetected++;
            }
            else {
                A[i - eatDetected] = A[i];
                B[i - eatDetected] = B[i];
            }
        }
        A[i - eatDetected] = A[i];
    } while (eatDetected);
    return (count + 1);
}
*/
/*
struct sNode
{
    char data;
    struct sNode* next;
};

void push(struct sNode** top_ref, int new_data) {
    struct sNode* new_node = (struct sNode*) malloc(sizeof(struct sNode));
    if (new_node == NULL)
        return;
    new_node->data = new_data;
    new_node->next = (*top_ref);
    (*top_ref) = new_node;
}

int pop(struct sNode** top_ref) {
    char res;
    struct sNode* top;
    if (*top_ref == NULL) {
        printf("Stack overflow n");
        getchar();
        exit(0);
    }
    else {
        top = *top_ref;
        res = top->data;
        *top_ref = top->next;
        free(top);
        return res;
    }
}

int isMatchingPair(char character1, char character2) {
    if (character1 == '(' && character2 == ')')
        return 1;
    else
        return 0;
}

int solution(char* S) {
    int i = 0;
    struct sNode* stack = NULL;
    while (S[i]) {
        if (S[i] == '(')
            push(&stack, S[i]);

        if (S[i] == ')') {
            if (stack == NULL)
                return 0;
            else if (!isMatchingPair(pop(&stack), S[i]))
                return 0;
        }
        i++;
    }

    if (stack == NULL)
        return 1;

    return 0;
}
*/
/*
int compare(const void* a, const void* b) {
    long long la = *(int*)a;
    long long lb = *(int*)b;
    if (la > lb)
        return 1;
    if (la < lb)
        return -1;
    return 0;
}

void myMemCpy(int* dest, int* src, int size)
{
    int* csrc = (int*)src;
    int* cdest = (int*)dest;

    for (int i = 0; i < size; i++)
        cdest[i] = csrc[i];
}

int AC[100000];
int solution(int A[], int N) {
    int val, count = 0;
    myMemCpy(AC, A, N);
    qsort(AC, N, sizeof(int), compare);
    val = AC[N / 2];
    
    for (int i = 0; i < N; i++) {
        if (val == AC[i]) {
            count++;
        }
    }

    if (count > (N / 2)) {
        for (int i = 0; i < N; i++) {
            if (val == A[i]) {
                return i;
            }
        }
    }
    return -1;
}
*/
/*
int compare(const void* a, const void* b) {
    long long la = *(int*)a;
    long long lb = *(int*)b;
    if (la > lb)
        return 1;
    if (la < lb)
        return -1;
    return 0;
}

void myMemCpy(int* dest, int* src, int size)
{
    int* csrc = (int*)src;
    int* cdest = (int*)dest;

    for (int i = 0; i < size; i++)
        cdest[i] = csrc[i];
}

int AC[100000];
int solution(int A[], int N) {
    int val, count = 0, j = 0, eq = 0;
    myMemCpy(AC, A, N);
    qsort(AC, N, sizeof(int), compare);
    val = AC[N / 2];

    for (int i = 0; i < N; i++) {
        if (val == AC[i]) {
            count++;
        }
    }

    for (int i = 0; i < (N - 1); i++) {
        if (val == A[i]) {
            j++;
        }
        if ((j > ((i + 1) / 2)) && ((count - j) > (N - i - 1) / 2)) {
            eq++;
        }
    }
    return eq;
}
*/
/*
int solution(int A[], int N) {
    int min = 200001;
    int max = 0;

    for (int i = 0; i < N; i++) {
        min = A[i] < min ? A[i] : min;
        max = A[i] - min > max ? A[i] - min : max;
    }
    return max;
}
*/

/*
int compare(const void* a, const void* b) {
    return (*(int*)a - *(int*)b);
}

int solution(int A[], int N) {

    if (N == 0)
        return 0;

    if (N == 1)
        return A[1];

    int sum = A[0] + A [1];

    qsort(A, N, sizeof(int), compare);
    for (int i = 2; i < N; i++) {
        sum += sum + A[i];
    }

    return sum;
}
*/
/*
char compressedstr[200000 + 1];
char newstr[200000 + 1];

int compressedStringLen(char* S) {
    int i = 0;
    int j = 0;
    int count = 0;

    memset(compressedstr, 0, strlen(S)+1);

    while (S[i] != 0)
    {
        if (S[i] == S[i + 1]) {
            count++;
        }
        else {
            if (count > 0) {
                count += 1;
                sprintf(&compressedstr[strlen(compressedstr)], "%d", count);
            }
            compressedstr[strlen(compressedstr)] = S[i];
            j++;
            count = 0;
        }
        i++;
    }

    return (int)strlen(compressedstr);
}

int solution(char* S, int K) {
    int origLen = compressedStringLen(S);
    int len;
    for(int i = K; i < (int)(strlen(S) - K); i++)
    {
        memcpy(newstr, S, i);
        memcpy(&newstr[i], &S[i + K], strlen(S) - i - K);
        len=compressedStringLen(newstr);
        if (len < origLen) {
            origLen = len;
        }
    }

    return origLen;
}
*/
/*
int solution(int A[], int N) {
    int sum = 0, tmpsum = 0;
    int lastVal = 0;
    int i = 0, neg = 0;

    sum = A[0];

    for (i = lastVal; i < N; i++) {
        if (sum > 0) {
            if (A[i] < 0) {
                if (i == lastVal) {
                    lastVal++;
                    neg = 0;
                    continue;
                }
                neg = 1;
            }
            else {
                if (i == lastVal) {
                    if (neg == 1) {
                        lastVal++;
                        continue;
                    }
                }
            }
        }

        tmpsum += A[i];
        if (sum < 0 && sum <= A[i]) {
            sum = A[i];
            lastVal = i;
            tmpsum = A[i];
        }
        else
        if (tmpsum > sum) {
            sum = tmpsum;
            if (lastVal != (N - 1) && i == (N - 1)) {
                if (neg == 0) {
                    break;
                }
                i = lastVal++;
                tmpsum = 0;
            }
        }
        else {
            if (i == (N - 1)) {
                i = lastVal++;
                tmpsum = 0;
            }
        }
    }

    if (sum < A[N - 1]) {
        sum = A[N - 1];
    }

    return sum;
}
*/

//char teststr[200000];
int main()
{
    int A[] = { -1, 2, -2, -4, 5, -4, 7};
//    int B[] = { 0, 1, 0, 0, 0 };
    //char A[] = "AAAABCCCDEEEEEEAAAAAAAAABCCCDEEEEEEAAAAAAAAABCCCDEEEEEEAAAAAAAAABCCCDEEEEEEAAAAAAAAABCCCDEEEEEEAAAAAAAAABCCCDEEEEEEAAAAAAAAABCCCDEEEEEEAAAAAAAAABCCCDEEEEEEAAAAAAAAABCCCDEEEEEEAAAAA";
    //for (int i = 0; i < 100; i++) {
    //    memcpy(&teststr[i * strlen(A)], A, strlen(A));
    //}
    
 //   std::cout << solution(A, sizeof(A)/sizeof(int)) <<"\n";
//    std::cout << solution(15) << "\n";
    
    return 0;
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu
