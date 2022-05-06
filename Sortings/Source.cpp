#include <cmath>
#include <cstdio>
#include <vector>
#include <iostream>
#include <algorithm>
#include <random>
#include <stack>
#include <queue>
#include <list>
#include <iostream>
#include <time.h>
#include <chrono> 
#include <assert.h> 
using namespace std;

int rnd() {
    return rand() % 20 - 10;
}

void printArr(int* arr, int size) {
    for (int i = 0; i < size; i++) {
        std::cout << arr[i] << " ";
    }
    std::cout << '\n';
}

void insertionSortDesc(int* arr, int n) {
    int value;
    for (int i = 1; i < n; i++) {
        value = arr[i];
        int j = i - 1;
        while (j >= 0 && value > arr[j]) {
            arr[j + 1] = arr[j];
            j--;
        }
        arr[j + 1] = value;
    }
}

void shellSort(int* arr, unsigned long size) {
    int gap = size / 2;
    int j, value;
    while (gap > 0) {
        for (unsigned long i = gap; i < size; i++) {
            value = arr[i];
            j = i;
            while (j >= gap && value < arr[j-gap]) {
                arr[j] = arr[j - gap];
                j -= gap;
            }
            arr[j] = value;
        }
        gap /= 2;
    }
}


/*
void shellSort(int* arr, int size) {
    int gap = size / 2;
    int j, value;
    while (gap > 0) {
        for (int i = gap; i < size; i++) {
            value = arr[i];
            j = i;
            while (j > 0 && arr[j - 1] < value) {
                arr[j] = arr[j - 1];
            }
            arr[j] = value;
        }
        gap /= 2;
    }
}

void shellSort(int * array, int n) {
    // Rearrange elements at each n/2, n/4, n/8, ... intervals
    for (int interval = n / 2; interval > 0; interval /= 2) {
        for (int i = interval; i < n; i += 1) {
            int temp = array[i];
            int j;
            for (j = i; j >= interval && array[j - interval] > temp; j -= interval) {
                array[j] = array[j - interval];
            }
            array[j] = temp;
        }
    }
}
*/



void shellSort2(int* arr, int size) {
    for (int gap = size / 2; gap > 0; gap/=2) {
        for (int i = gap; i < size; i++) {
            int value = arr[i];
            int j = i - gap;
            while (value < arr[j] && j >= 0) {
                arr[j + gap] = arr[j];
                j -= gap;
            }
            arr[j + gap] = value;
        }
    }
}


void binarySearch(int* arr, int n, int element) {
    int start = 0;
    int end = n - 1;
    int index;
    while (start <= end) {
        index = (start + end) / 2;
        if (arr[index] == element) {
            std::cout << "Element is on index " << index << ".\n";
            return;
        }
        else {
            if (arr[index] > element) {
                end = index - 1;
            }
            else {
                start = index + 1;
            }
        }
    }
    std::cout << "The element was not found.\n";
}

void binarySearchAscending(int* arr, int n, int element) {
    int start = 0;
    int end = n - 1;
    int index;
    while (start <= end) {
        index = (start + end) / 2;
        if (arr[index] == element) {
            std::cout << "Element is on index " << index << ".\n";
            return;
        }
        else {
            if (arr[index] > element) {
                end = index - 1;
            }
            else {
                start = index + 1;
            }
        }
    }
    std::cout << "The element was not found.\n";
}

void binarySearchDescending(int* arr, int n, int element) {
    int start = 0;
    int end = n - 1;
    int index;
    while (start <= end) {
        index = (start + end) / 2;
        if (arr[index] == element) {
            std::cout << "Element is on index " << index << ".\n";
            return;
        }
        else {
            if (arr[index] > element) {
                end = index - 1;
            }
            else {
                start = index + 1;
            }
        }
    }
    std::cout << "The element was not found.\n";
}

int binarySearchIndex(int* arr, int size,
    int element) {
    int start = 0;
    int end = size - 1;
    int index = (start + end) / 2;
    if (element != arr[index]) {

        while (start < end) {
            if (element > arr[index]) {
                start = index + 1;
            }
            else {
                end = index - 1;
            }
            index = (start + end) / 2;
        }

        if (element > arr[index]) {
            index++;
        }
    }
    return index;
}

void binarySearchRecursive(int* arr, int n, int element, int start, int end) {
    int index = (start + end) / 2;
    if (arr[index] == element) {
        std::cout << "Element is on index " << index << ".\n";
        return;
    }
    if (start >= end) {
        std::cout << "The element was not found.\n";
        return;
    }

    if (arr[index] > element) {
        binarySearchRecursive(arr, n, element, start, index - 1);
    }
    else {
        binarySearchRecursive(arr, n, element, index + 1, end);
    }

}

void quickSort(int* arr, int start, int end) {
    if (start >= end) {
        return;
    }

    std::random_device                  rand_dev;
    std::mt19937                        generator(rand_dev());
    std::uniform_int_distribution<int>  distr(start, end);
    int pivot = arr[distr(generator)];
    int leftIterator = start;
    int rightIterator = end;

    while (leftIterator <= rightIterator) {
        while (arr[leftIterator] < pivot) {
            ++leftIterator;
        }
        while (arr[rightIterator] > pivot) {
            --rightIterator;
        }
        if (leftIterator <= rightIterator) {
            std::swap(arr[leftIterator], arr[rightIterator]);
            ++leftIterator;
            --rightIterator;
        }
    }
    quickSort(arr, start, rightIterator);
    quickSort(arr, leftIterator, end);
}

void quickSortIteratively(int* arr, int n) {
    std::queue<int> boundaryHolder;
    int start, end, leftIterator, rightIterator, pivot;
    boundaryHolder.push(0);
    boundaryHolder.push(n - 1);
    while (!boundaryHolder.empty()) {
        leftIterator = start = boundaryHolder.front();
        boundaryHolder.pop();
        rightIterator = end = boundaryHolder.front();
        boundaryHolder.pop();
        pivot = arr[(leftIterator + rightIterator) / 2];
        while (leftIterator <= rightIterator) {
            while (arr[leftIterator] < pivot) {
                ++leftIterator;
            }
            while (arr[rightIterator] > pivot) {
                --rightIterator;
            }
            if (leftIterator <= rightIterator) {
                std::swap(arr[leftIterator], arr[rightIterator]);
                ++leftIterator;
                --rightIterator;
            }
        }
        if (start < rightIterator) {
            boundaryHolder.push(start);
            boundaryHolder.push(rightIterator);
        }
        if (leftIterator < end) {
            boundaryHolder.push(leftIterator);
            boundaryHolder.push(end);
        }
    }
}









//MERGE FUNCTION
void merge(int* arr, int* bfarr, int start, int index, int end) {
    int leftArrayIterator = start;
    int rightArrayIterator = index;
    int bfIterator = start;
    while (leftArrayIterator < index && rightArrayIterator <= end) {
        if (arr[leftArrayIterator] < arr[rightArrayIterator]) {
            bfarr[bfIterator] = arr[leftArrayIterator];
            ++leftArrayIterator;
            ++bfIterator;
        }
        else {
            bfarr[bfIterator] = arr[rightArrayIterator];
            ++rightArrayIterator;
            ++bfIterator;
        }
    }
    for (int k = leftArrayIterator; k < index; k++) {
        bfarr[bfIterator] = arr[k];
        ++bfIterator;
    }
    for (int k = rightArrayIterator; k <= end; k++) {
        bfarr[bfIterator] = arr[k];
        ++bfIterator;
    }
    for (int k = start; k <= end; k++) {
        arr[k] = bfarr[k];
    }
}
//////////////////







void partition(int* arr, int* bfarr, int start, int end) {
    int index = (start + end) / 2;
    if (start >= end) {
        return;
    }
    partition(arr, bfarr, start, index);
    partition(arr, bfarr, index + 1, end);
    merge(arr, bfarr, start, index + 1, end);
}

void MergeSort(int* arr, int size) {
    int* bfarr = new int[size];
    partition(arr, bfarr, 0, size - 1);
    delete[] bfarr;
}

void MergeSort2(int* arr, int start, int end) {
    int* bfarr = new int[end-start+1];
    partition(arr, bfarr, start, end);
    delete[] bfarr;
}

void mergeSortIteratively(int* arr, int size) {
    int* bfarr = new int[size];
    for (int i = 1; i < size; i *= 2) {
        for (int j = 0; j < size; j += i * 2) {
            merge(arr, bfarr, j, min(j + i, size - 1), min(j + i * 2 - 1, size - 1));
        }
    }
}

void mergeSortIterativelyVector(int* arr, int size) {
    int* bfarr = new int[size];

    std::vector<std::vector<int>> partitions;
    partitions.push_back(std::vector<int> {0, size - 1});
    for (int i = 0; i < partitions.size(); i++) {
        int index = (partitions[i][0] + partitions[i][1]) / 2;
        if (partitions[i][0] < index) {
            partitions.push_back(std::vector<int> {partitions[i][0], index});
        }
        if (index + 1 < partitions[i][1]) {
            partitions.push_back(std::vector<int> {index + 1, partitions[i][1]});
        }
    }

    while (!partitions.empty()) {
        int index = (partitions.back()[0] + partitions.back()[1]) / 2;
        merge(arr, bfarr, partitions.back()[0], index + 1, partitions.back()[1]);
        partitions.pop_back();
    }
}

void mergeSortIterativelyList(int* arr, int size) {
    int* bfarr = new int[size];
    int index;

    std::list <std::pair<int, int>> partitions;
    partitions.push_back(std::pair<int, int>(0, size - 1));

    std::list<std::pair<int, int>>::iterator it;

    for (it = partitions.begin(); it != partitions.end(); ++it) {
        index = ((*it).first + (*it).second) / 2;
        if ((*it).first < index) {
            partitions.push_back(std::pair<int, int>((*it).first, index));
        }
        if (index + 1 < (*it).second) {
            partitions.push_back(std::pair<int, int>(index + 1, (*it).second));
        }
    }

    while (!partitions.empty()) {
        index = (partitions.back().first + partitions.back().second) / 2;
        merge(arr, bfarr, partitions.back().first, index + 1, partitions.back().second);
        partitions.pop_back();
    }
}






//HEAP
void maxHeapify(int* arr, int size, int rootIndex)
{
    int index = rootIndex;
    int left = 2 * index + 1;
    int right = 2 * index + 2;

    if (left<size && arr[left]>arr[index]) {
        index = left;
    }

    if (right<size && arr[right]>arr[index]) {
        index = right;
    }

    if (index != rootIndex) {
        std::swap(arr[rootIndex], arr[index]);
        maxHeapify(arr, size, index);
    }
}

void buildHeap(int* arr, int size)
{
    int start = (size / 2) - 1;
    for (int i = start; i >= 0; i--) {
        maxHeapify(arr, size, i);
        //minHeapify(arr, size, i);
    }
}


/////////////////////////////////////////////////////






//////////////////////////////INTRO SORT/////////////////////
/*
////////////////////////////////////////////////////////////////////
/*
template <class _RanIt, class _Pr>
_CONSTEXPR20 void _Make_heap_unchecked(_RanIt _First, _RanIt _Last, _Pr _Pred) {
    // make [_First, _Last) into a heap
    using _Diff   = _Iter_diff_t<_RanIt>;
    _Diff _Bottom = _Last - _First;
    for (_Diff _Hole = _Bottom >> 1; _Hole > 0;) { // shift for codegen
        // reheap top half, bottom to top
        --_Hole;
        _Iter_value_t<_RanIt> _Val = _STD move(*(_First + _Hole));
        _Pop_heap_hole_by_index(_First, _Hole, _Bottom, _STD move(_Val), _Pred);
    }
}
*/

void buildHeapIterativelyMax(int* arr, int size) {
    int start = size / 2 - 1, index, lastIndex, left, right;
    for (int i = start; i >= 0; i--) {
        lastIndex = i;
        while (true) {
            index = lastIndex;
            left = 2 * index + 1;
            right = 2 * index + 2;

            if (left<size && arr[left]>arr[index]) {
                index = left;
            }

            if (right<size && arr[right]>arr[index]) {
                index = right;
            }

            if (index == lastIndex) {
                break;
            }
            std::swap(arr[lastIndex], arr[index]);
            lastIndex = index;
        }
    }
}


void siftDownIterative(int* arr, int start, int end) {
    int index, left, right;
    while (true) {
        index = start;
        left = 2 * index + 1;
        if (left <= end) {
            if (arr[left] > arr[index]) {
                index = left;
            }
            right = left + 1;
            if (right <= end && arr[right] > arr[index]) {
                index = right;
            }
        }
        if (index == start) {
            break;
        }
        std::swap(arr[start], arr[index]);
        start = index;
    }
}



void heapSortIterative(int* arr, int size) {
    buildHeapIterativelyMax(arr, size);
    int end = size-1;
    while (end > 0) {
        std::swap(arr[0], arr[end]);
        end -= 1;
        siftDownIterative(arr, 0, end);
    }
}


//////////////////////////INTRO SORT ITERATIVE///////////////////////////
/*
*/
//////////////////////////INTRO SORT ITERATIVE///////////////////////////

/*
void compareQuickRecIntroRec(int* arr, int size) {
    int* arrcpy1 = new int[size];
    int* arrcpy2 = new int[size];

    for (int i = 0; i < size; i++) {
        arrcpy1[i] = arr[i];
        arrcpy2[i] = arr[i];
    }


    //RECURSIVE QUICKSORT arrcpy1
    std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();
    quickSort(arrcpy1, 0, size - 1);
    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
    std::cout << "QuickSortRecursive:\n Elapsed nanoseconds:"
        << std::chrono::duration_cast<std::chrono::nanoseconds> (t1 - t0).count()
        << std::endl;

    //RECURSIVE INTROSORT arrcpy3
    //RECURSION doesn't work for 2000000 elements
    //Dang 2000000 elements are a lot of elements
    std::chrono::steady_clock::time_point t4 = std::chrono::steady_clock::now();
    introSort(arrcpy2, size);
    std::chrono::steady_clock::time_point t5 = std::chrono::steady_clock::now();
    std::cout << "IntroSortRecursive:\n Elapsed nanoseconds:"
        << std::chrono::duration_cast<std::chrono::nanoseconds> (t5 - t4).count()
        << std::endl;

    bool equal = true;
    for (int i = 0; i < size; i++) {
        std::cout << arrcpy2[i] << " ";
        if (arrcpy1[i] != arrcpy2[i]) {
            equal = false;
        }
    }
    std::cout << "\nArrays qIt and introIt are equal: " << equal << "\n";
}
*/

/*
void compareHeapQuick(int* arr, int size) {
    int* arrcpy1 = new int[size];
    int* arrcpy2 = new int[size];
    int* arrcpy3 = new int[size];

    for (int i = 0; i < size; i++) {
        arrcpy1[i] = arr[i];
        arrcpy2[i] = arr[i];
        arrcpy3[i] = arr[i];
    }

    quickSort(arrcpy1, 40, size - 21);
    heapSortIntro(arrcpy2, 40, size - 21);
    //quickSortIteratively(arrcpy3, size);

    bool equal = true;
    //bool equal2 = true;
    for (int i = 0; i < size; i++) {
        if (arrcpy1[i] != arrcpy2[i]) {
            equal = false;
        }
        //if (arrcpy1[i] != arrcpy3[i]) {
        //    equal2 = false;
        //}
    }
    std::cout << "Quick equal to Heap Intro: " << equal<< "\n";
    //std::cout << "Quick equal to Quick Iterative: " << equal2 << "\n";

    assert(equal);
    //assert(equal2);
}
*/



////////////////////////////INSERTION SORT/////////////////////////////////
void insertionSortIntro(int* arr, int start, int end) {
    int value;
    int j;
    for (int i = start + 1; i <= end; i++) {
        value = arr[i];
        j = i;
        while (value<arr[j-1]&&j>start) {
            arr[j] = arr[j-1];
            j--;
        }
        arr[j] = value;
    }
}


///////////////////////SHELL SORT/////////////////////////
void shellSortIntro(int* arr, int start, int end) {
    int gap = (end-start)/2;
    while (gap > 0) {
        int j, value;
        int i = start+gap;
        while (i <= end) {
            j = i;
            value = arr[i];
            while (value < arr[j-gap] && j >= start+gap) {
                arr[j] = arr[j-gap];
                j -= gap;
            }
            arr[j] = value;
            ++i;
        }
        gap /= 2;
    }
}

///////////////////////////////////////////////////////////////////////////



///////////////////////////MODIFIED HEAP AND HEAPSORT//////////////////////////////////
//MODIFIED HEAP
inline void maxHeapifyIntro(int* arr, int startingIndex, int start, int end, int heapDistance) {
    int index = startingIndex;
    int left = start + (heapDistance * 2) + 1;
    int right = start + (heapDistance * 2) + 2;

    if (left <= end && arr[index] < arr[left]) {
        heapDistance = left - start;
        index = left;
    }

    if (right <= end && arr[index] < arr[right]) {
        heapDistance = right - start;
        index = right;
    }
    if (index != startingIndex) {
        std::swap(arr[index], arr[startingIndex]);
        maxHeapifyIntro(arr, index, start, end, heapDistance);
    }
}


inline void buildHeapIntro(int* arr, int start, int end) {
    int size = end - start + 1;
    int distance = (size / 2) - 1;
    for (; distance >= 0; distance--) {
        maxHeapifyIntro(arr, start + distance, start, end, distance);
    }
}

inline void siftDownIntro(int* arr, int startingIndex, int start, int end, int heapDistance) {
    int index = startingIndex;
    int left = start + (heapDistance * 2) + 1;
    int right = start + (heapDistance * 2) + 2;

    if (left <= end && arr[index] < arr[left]) {
        heapDistance = left - start;
        index = left;
    }

    if (right <= end && arr[index] < arr[right]) {
        index = right;
        heapDistance = right - start;
    }
    if (index != startingIndex) {
        std::swap(arr[index], arr[startingIndex]);
        siftDownIntro(arr, index, start, end, heapDistance);
    }
}

inline void heapSortIntro(int* arr, int start, int end) {
    int endIterator = end;
    buildHeapIntro(arr, start, end);
    while (endIterator > start) {
        std::swap(arr[start], arr[endIterator]);
        --endIterator;
        siftDownIntro(arr, start, start, endIterator, 0);
    }
}
///////////////////////////////HEAP SORT INTRO//////////////////////////////


///////////////////////////////INTRO SORT RECURSIVE//////////////////////////////
inline void introSortProcedure(int* arr, int start, int end, int  depth) {
    //if (start >= end) {
    //    return;
    //}

    if (end - start < 144) {
        insertionSortIntro(arr, start, end);
    }

    else if (depth == 0) {
        heapSortIntro(arr, start, end);
    }

    else {
        int pivot = arr[(start + end) / 2];
        int leftIterator = start;
        int rightIterator = end;
        while (leftIterator <= rightIterator) {
            while (arr[leftIterator] < pivot) {
                ++leftIterator;
            }
            while (arr[rightIterator] > pivot) {
                --rightIterator;
            }
            if (leftIterator <= rightIterator) {
                std::swap(arr[leftIterator], arr[rightIterator]);
                ++leftIterator;
                --rightIterator;
            }
        }
        introSortProcedure(arr, start, rightIterator, depth - 1);
        introSortProcedure(arr, leftIterator, end, depth - 1);
    }

}

inline void introSort(int* arr, int size) {
    int depth = log2(size) * 2;
    introSortProcedure(arr, 0, size - 1, depth);
}



void introSortProcedure2(int* arr, int start, int end, int  depth) {
    if (start >= end) {
        return;
    }

    if (end - start < 900) {
        if (end - start < 100) {
            insertionSortIntro(arr, start, end);
        }
        else shellSortIntro(arr, start, end);
    }
    else if (depth == 0) {
        heapSortIntro(arr, start, end);
    }

    else {
        int pivot = arr[(start + end) / 2];
        int leftIterator = start;
        int rightIterator = end;
        while (leftIterator <= rightIterator) {
            while (arr[leftIterator] < pivot) {
                ++leftIterator;
            }
            while (arr[rightIterator] > pivot) {
                --rightIterator;
            }
            if (leftIterator <= rightIterator) {
                std::swap(arr[leftIterator], arr[rightIterator]);
                ++leftIterator;
                --rightIterator;
            }
        }
        introSortProcedure2(arr, start, rightIterator, depth - 1);
        introSortProcedure2(arr, leftIterator, end, depth - 1);
    }

}

void introSort2(int* arr, int size) {
    int depth = log(size) * 2;
    introSortProcedure2(arr, 0, size - 1, depth);
}


///////////////////////////ITERATIVE VERSIONS////////////////////////////////////


/*
void heapSortProcedureIntroIterative(int* arr, int start, int end) {
    int endIterator = end;
    int index, left, right, heapDistance, startingIndex;
    while (endIterator > start) {
        std::swap(arr[start], arr[endIterator]);
        --endIterator;
        startingIndex = start;
        heapDistance = 0;
        while (true) {
            index = startingIndex;
            left = start + heapDistance * 2 + 1;
            right = start + heapDistance * 2 + 2;

            if (left <= end && arr[index] < arr[left]) {
                heapDistance = left - start;
                index = left;
            }

            if (right <= end && arr[index] < arr[right]) {
                heapDistance = right - start;
                index = right;
            }

            if (index == startingIndex) {
                break;
            }

            std::swap(arr[index], arr[startingIndex]);
            startingIndex = index;
        }
    }
}*/

void buildHeapMaxIntroIterative(int* arr, int start, int end) {
    int size = end - start + 1;
    int distance = size / 2 - 1;
    int index, lastIndex, left, right, heapDistance;
    for (; distance >= 0; distance--) {
        lastIndex = start + distance;
        heapDistance = distance;
        while (true) {
            index = lastIndex;
            left = start + heapDistance * 2 + 1;
            right = start + heapDistance * 2 + 2;

            if (left <= end && arr[index] < arr[left]) {
                heapDistance = left-start;
                index = left;
            }

            if (right <= end && arr[index] < arr[right]) {
                heapDistance = right - start;
                index = right;
            }

            if (index == lastIndex) {
                break;
            }

            std::swap(arr[index], arr[lastIndex]);
            lastIndex = index;
        }
    }
}


void siftDownIntroIterative(int* arr, int start, int end) {
    int index, left, right;
    int heapDistance = 0;
    int startingIndex = start;
    while (true) {
        index = startingIndex;
        left = start + heapDistance * 2 + 1;
        right = start + heapDistance * 2 + 2;

        if (left <= end && arr[index] < arr[left]) {
            heapDistance = left - start;
            index = left;
        }

        if (right <= end && arr[index] < arr[right]) {
            heapDistance = right - start;
            index = right;
        }

        if (index == startingIndex) {
            break;
        }
        std::swap(arr[index], arr[startingIndex]);
        startingIndex = index;
    }
}

void heapSortIntroIterative(int* arr, int start, int end) {
    buildHeapMaxIntroIterative(arr, start, end);
    int endIterator = end;
    while (endIterator > start) {
        std::swap(arr[start], arr[endIterator]);
        --endIterator;
        siftDownIntroIterative(arr, start, endIterator);
    }
}

/*
void introSortIterativeQueue(int* arr, int size) {
    //std::list <std::queue<int>> partitions;
    std::queue<std::queue<int>> partitions;
    std::queue<int> initialQ;
    initialQ.push(0);
    initialQ.push(size - 1);
    partitions.push(initialQ);
    //std::list<std::vector<int>>::iterator it; 
    int depth = log(size) * 2;
    int leftIterator, start, end, rightIterator, pivot;
    while (!partitions.empty()) {
        std::queue<int> newBoundaries;
        while (!partitions.front().empty()) {
            leftIterator = start = partitions.front().front();
            partitions.front().pop();
            rightIterator = end = partitions.front().front();
            partitions.front().pop();

            if (end - start < 64) {
                insertionSortIntro(arr, start, end);
            }
            else if (depth == 0) {
                heapSortIntroIterative(arr, start, end);
            }
            else {
                pivot = arr[(start + end) / 2];
                while (leftIterator <= rightIterator) {
                    while (arr[leftIterator] < pivot) {
                        ++leftIterator;
                    }
                    while (arr[rightIterator] > pivot) {
                        --rightIterator;
                    }
                    if (leftIterator <= rightIterator) {
                        std::swap(arr[leftIterator], arr[rightIterator]);
                        ++leftIterator;
                        --rightIterator;
                    }
                }

                if (start < rightIterator) {
                    newBoundaries.push(start);
                    newBoundaries.push(rightIterator);
                }

                if (leftIterator < end) {
                    newBoundaries.push(leftIterator);
                    newBoundaries.push(end);
                }
            }
        }
        partitions.pop();
        if (!newBoundaries.empty()) {
            --depth;
            partitions.push(newBoundaries);
        }
    }
}
*/


//CURRENT INTROSORT ITERATIVE QUEUE!!!
/*
void introSortIterativeQueue(int* arr, int size) {
    //std::list <std::queue<int>> partitions;
    std::queue<int> partitions;
    std::queue<int> bufferPartitions;
    partitions.push(0);
    partitions.push(size - 1);
    //std::list<std::vector<int>>::iterator it; 
    int depth = log(size)*2;
    int leftIterator, start, end, rightIterator, pivot;
    while (!partitions.empty()||!bufferPartitions.empty()) {
        while (!partitions.empty()) {
            leftIterator = start = partitions.front();
            partitions.pop();
            rightIterator = end = partitions.front();
            partitions.pop();

            if (end - start < 100) {
                insertionSortIntro(arr, start, end);
            }
            else if (depth < 1) {
                //std::cout << "heapSort\n";
                heapSortIntroIterative(arr, start, end);
            }
            else {
                pivot = arr[(start + end) / 2];
                while (leftIterator <= rightIterator) {
                    while (arr[leftIterator] < pivot) {
                        ++leftIterator;
                    }
                    while (arr[rightIterator] > pivot) {
                        --rightIterator;
                    }
                    if (leftIterator <= rightIterator) {
                        std::swap(arr[leftIterator], arr[rightIterator]);
                        ++leftIterator;
                        --rightIterator;
                    }
                }

                if (start < rightIterator) {
                    bufferPartitions.push(start);
                    bufferPartitions.push(rightIterator);
                }

                if (leftIterator < end) {
                    bufferPartitions.push(leftIterator);
                    bufferPartitions.push(end);
                }
            }
        }
        --depth;
        while (!bufferPartitions.empty()) {
            leftIterator = start = bufferPartitions.front();
            bufferPartitions.pop();
            rightIterator = end = bufferPartitions.front();
            bufferPartitions.pop();

            if (end - start < 100) {
                insertionSortIntro(arr, start, end);
            }
            else if (depth < 1) {
                //std::cout << "heapSort\n";
                heapSortIntroIterative(arr, start, end);
            }
            else {
                pivot = arr[(start + end) / 2];
                while (leftIterator <= rightIterator) {
                    while (arr[leftIterator] < pivot) {
                        ++leftIterator;
                    }
                    while (arr[rightIterator] > pivot) {
                        --rightIterator;
                    }
                    if (leftIterator <= rightIterator) {
                        std::swap(arr[leftIterator], arr[rightIterator]);
                        ++leftIterator;
                        --rightIterator;
                    }
                }

                if (start < rightIterator) {
                    partitions.push(start);
                    partitions.push(rightIterator);
                }

                if (leftIterator < end) {
                    partitions.push(leftIterator);
                    partitions.push(end);
                }
            }
        }
        --depth;
    }
}
*/

inline void introSortIterativeQueue(int* arr, int size) {
    if (size > 1) {
        std::queue<int> bufferQueue1, bufferQueue2, * partitions = &bufferQueue1, * nextPartition = &bufferQueue2, * tmp;
        partitions->push(0);
        partitions->push(size - 1);
        int depth = log(size) * 2;
        int leftIterator, start, end, rightIterator, pivot; 
        while (!partitions->empty()) {
            while (!partitions->empty()) {
                leftIterator = start = partitions->front();
                partitions->pop();
                rightIterator = end = partitions->front();
                partitions->pop();
                if (end - start < 144) {
                    insertionSortIntro(arr, start, end);
                }
                else if (depth == 0) {
                    heapSortIntroIterative(arr, start, end);
                }
                else {
                    pivot = arr[(start + end)/2];
                    while (leftIterator <= rightIterator) {
                        while (arr[leftIterator] < pivot) {
                            ++leftIterator;
                        }
                        while (arr[rightIterator] > pivot) {
                            --rightIterator;
                        }
                        if (leftIterator <= rightIterator) {
                            std::swap(arr[leftIterator], arr[rightIterator]);
                            ++leftIterator;
                            --rightIterator;
                        }
                    }
                    if (start < rightIterator) {
                        nextPartition->push(start);
                        nextPartition->push(rightIterator);
                    }
                    if (leftIterator < end) {
                        nextPartition->push(leftIterator);
                        nextPartition->push(end);
                    }
                }
            }
            std::swap(partitions, nextPartition);
            --depth;
        }
    }
}


inline void introSortIterativeQueueLong(int* arr, int size) {
    if (size > 1) {
        std::queue<int> bufferQueue1, bufferQueue2, * partitions = &bufferQueue1, * nextPartition = &bufferQueue2, * tmp;
        partitions->push(0);
        partitions->push(size - 1);
        int depth = log2(size)*2;
        int leftIterator, start, end, rightIterator, pivot;
        int subSize, heapDistance, distance, index;
        while (!partitions->empty()) {
            while (!partitions->empty()) {
                leftIterator = start = partitions->front();
                partitions->pop();
                rightIterator = end = partitions->front();
                partitions->pop();
                subSize = end - start;
                if (subSize < 144) {
                    insertionSortIntro(arr, start, end);
                }
                else if (depth == 0) { 
                    distance = subSize / 2 - 1;
                    heapDistance = distance;
                    for (; distance >= 0; distance--) {
                        pivot = start + distance;
                        while (true) {
                            index = pivot;
                            leftIterator = start + heapDistance * 2 + 1;
                            rightIterator = start + heapDistance * 2 + 2;

                            if (leftIterator <= end && arr[index] < arr[leftIterator]) {
                                heapDistance = leftIterator - start;
                                index = leftIterator;
                            }

                            if (rightIterator <= end && arr[index] < arr[rightIterator]) {
                                heapDistance = rightIterator - start;
                                index = rightIterator;
                            }

                            if (index == pivot) {
                                break;
                            }

                            std::swap(arr[index], arr[pivot]);
                            pivot = index;
                            --heapDistance;
                        }
                    }
                    while (end > start) {
                        std::swap(arr[start], arr[end]);
                        --end;
                        heapDistance = 0;
                        pivot = start;
                        while (true) {
                            index = pivot;
                            leftIterator = start + heapDistance * 2 + 1;
                            rightIterator = start + heapDistance * 2 + 2;

                            if (leftIterator <= end && arr[index] < arr[leftIterator]) {
                                heapDistance = leftIterator - start;
                                index = leftIterator;
                            }

                            if (rightIterator <= end && arr[index] < arr[rightIterator]) {
                                heapDistance = rightIterator - start;
                                index = rightIterator;
                            }

                            if (index == pivot) {
                                break;
                            }
                            std::swap(arr[index], arr[pivot]);
                            pivot = index;
                        }
                    }
                }
                else {
                    pivot = arr[(start + end) / 2];
                    while (leftIterator <= rightIterator) {
                        while (arr[leftIterator] < pivot) {
                            ++leftIterator;
                        }
                        while (arr[rightIterator] > pivot) {
                            --rightIterator;
                        }
                        if (leftIterator <= rightIterator) {
                            std::swap(arr[leftIterator], arr[rightIterator]);
                            ++leftIterator;
                            --rightIterator;
                        }
                    }
                    if (start < rightIterator) {
                        nextPartition->push(start);
                        nextPartition->push(rightIterator);
                    }
                    if (leftIterator < end) {
                        nextPartition->push(leftIterator);
                        nextPartition->push(end);
                    }
                }
            }
            std::swap(partitions, nextPartition);
            --depth;
        }
    }
}


void introSortIterativeList(int* arr, int size) {
    std::list<std::vector<int>> partitions;
    partitions.push_back(std::vector<int>{0,size-1});
    std::list<std::vector<int>>::iterator it;
    int leftIterator, rightIterator, pivot;
    int depth = log(size) * 2;

    for (it = partitions.begin(); it != partitions.end(); ++it) {
        std::vector<int> boundaries;
        for (int i = 0; i < (*it).size(); i+=2) {
            if ((*it)[i + 1] - (*it)[i] < 16) {
                insertionSortIntro(arr, (*it)[i], (*it)[i + 1]);
            }
            else if (depth == 0) {
                heapSortIntroIterative(arr, (*it)[i], (*it)[i + 1]);
            }
            else {
                leftIterator = (*it)[i];
                rightIterator = (*it)[i + 1];
                //std::cout << "left = " << leftIterator << " right = " << rightIterator << "\n";
                pivot = arr[(leftIterator + rightIterator) / 2];
                while (leftIterator <= rightIterator) {
                    while (arr[leftIterator] < pivot) {
                        ++leftIterator;
                    }
                    while (arr[rightIterator] > pivot) {
                        --rightIterator;
                    }
                    if (leftIterator <= rightIterator) {
                        std::swap(arr[leftIterator], arr[rightIterator]);
                        ++leftIterator;
                        --rightIterator;
                    }
                }

                if ((*it)[i] < rightIterator) {
                    boundaries.push_back((*it)[i]);
                    boundaries.push_back(rightIterator);
                }

                if (leftIterator < (*it)[i+1]) {
                    boundaries.push_back(leftIterator);
                    boundaries.push_back((*it)[i + 1]);
                }
            }
        }
        if (!boundaries.empty()) {
            partitions.push_back(boundaries);
            --depth;
        }
    }
}


void introSortIterative(int* arr, int size) {
    std::queue<int> boundaryHolder;
    //std::queue<int> depth;
    int depth = (int)log(size)<<4;
    //std::cout << depth << "\n";
    //std::cout << depth << "\n";
    boundaryHolder.push(0);
    boundaryHolder.push(size - 1);
    //depth.push(log(size) * 2);
    int start, end, leftIterator, rightIterator, pivot;
    while (!boundaryHolder.empty()) {
        //std::cout << depth.front()<<"\n";
        leftIterator = start = boundaryHolder.front();
        boundaryHolder.pop();
        rightIterator = end = boundaryHolder.front();
        boundaryHolder.pop(); 
        //std::cout << depth << "\n";

        if (end - start <= 32) {
            //std::cout << "descending into insertion sort\n";
            insertionSortIntro(arr, start, end);
        }

        else if (depth == 0) {
            //std::cout << "descending into heap sort\n";
            heapSortIntro(arr, start, end);
        }
        else {
            //std::cout<<"descending into quick sort:\n";
            pivot = arr[(leftIterator + rightIterator) / 2];
            while (leftIterator <= rightIterator) {
                while (arr[leftIterator] < pivot) {
                    ++leftIterator;
                }
                while (arr[rightIterator] > pivot) {
                    --rightIterator;
                }
                if (leftIterator <= rightIterator) {
                    std::swap(arr[leftIterator], arr[rightIterator]);
                    ++leftIterator;
                    --rightIterator;
                }
            }

            if (start < rightIterator) {
                boundaryHolder.push(start);
                boundaryHolder.push(rightIterator);
                //depth.push(depth.front() - 1);
                //++depth;
            }

            if (leftIterator < end) {
                boundaryHolder.push(leftIterator);
                boundaryHolder.push(end);
                //depth.push(depth.front() - 1);
                //++depth;
            }
            --depth;
        }
    }
}


void compareQuickToIntro(int* arr, int size) {
    int* arrcpy1 = new int[size];
    int* arrcpy2 = new int[size];
    int* arrcpy3 = new int[size];
    std::vector<int> arrcpy4(size);
    int* arrcpy5 = new int[size];
    int* arrcpy6 = new int[size];

    for (int i = 0; i < size; i++) {
        arrcpy1[i] = arr[i];
        arrcpy2[i] = arr[i];
        arrcpy3[i] = arr[i];
        arrcpy4[i] = arr[i];
        arrcpy5[i] = arr[i];
        arrcpy6[i] = arr[i];
    }



    //RECURSIVE QUICKSORT arrcpy1
    std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();

    quickSort(arrcpy1, 0, size - 1);

    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
    std::cout << "QuickSortRecursive:\n Elapsed nanoseconds:"
        << std::chrono::duration_cast<std::chrono::nanoseconds> (t1 - t0).count()
        << std::endl;



    //ITERATIVE QUICKSORT arrcpy2
    std::chrono::steady_clock::time_point t6 = std::chrono::steady_clock::now();

    quickSortIteratively(arrcpy2, size);


    std::chrono::steady_clock::time_point t7 = std::chrono::steady_clock::now();
    std::cout << "QuickSortIterative:\n Elapsed nanoseconds:"
        << std::chrono::duration_cast<std::chrono::nanoseconds> (t7 - t6).count()
        << std::endl;



    //RECURSIVE INTROSORT arrcpy3
    //RECURSION doesn't work for 2000000 elements
    //Dang 2000000 elements are a lot of elements
    std::chrono::steady_clock::time_point t4 = std::chrono::steady_clock::now();

    introSort(arrcpy3, size);

    std::chrono::steady_clock::time_point t5 = std::chrono::steady_clock::now();
    std::cout << "My Sort:\n Elapsed nanoseconds:"
        << std::chrono::duration_cast<std::chrono::nanoseconds> (t5 - t4).count()
        << std::endl;

    bool equal = true;
    for (int i = 0; i < size; i++) {
        if (arrcpy2[i] != arrcpy3[i]) {
            equal = false;
            break;
        }
    }

    std::cout << "My Array Equal to Quick Array: " << equal << "\n";
    assert(equal);

    std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
    std::sort(arrcpy4.begin(), arrcpy4.end());
    std::chrono::steady_clock::time_point t3 = std::chrono::steady_clock::now();
    std::cout << "STL Sort:\n Elapsed nanoseconds:"
        << std::chrono::duration_cast<std::chrono::nanoseconds> (t3 - t2).count()
        << std::endl;

    bool equal2 = true;
    for (int i = 0; i < size; i++) {
        if (arrcpy4[i] != arrcpy3[i]) {
            equal2 = false;
            break;
        }
    }

    std::cout << "Stl Sort equal to My Sort1: " << equal2 << "\n";
    assert(equal2);

    //ITERATIVE INTROSORT arrcpy6
    std::chrono::steady_clock::time_point t12 = std::chrono::steady_clock::now();

    introSortIterativeQueue(arrcpy6, size);

    std::chrono::steady_clock::time_point t13 = std::chrono::steady_clock::now();
    std::cout << "Intro Iterative Queue:\n Elapsed nanoseconds:"
        << std::chrono::duration_cast<std::chrono::nanoseconds> (t13 - t12).count()
        << std::endl;

    bool equal4 = true;
    for (int i = 0; i < size; i++) {
        if (arrcpy4[i] != arrcpy6[i]) {
            equal4 = false;
            break;
        }
    }
        std::cout << "Stl Sort equal to My Sort2: " << equal4 << "\n";
        assert(equal4);
    /*
    //ITERATIVE INTROSORT arrcpy4
    std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
    introSortIterative(arrcpy4, size);
    std::chrono::steady_clock::time_point t3 = std::chrono::steady_clock::now();
    std::cout << "IntroSortIterative:\n Elapsed nanoseconds:"
        << std::chrono::duration_cast<std::chrono::nanoseconds> (t3 - t2).count()
        << std::endl;


    bool equal = true;
    for (int i = 0; i < size; i++) {
        if (arrcpy1[i] != arrcpy4[i]) {
            equal = false;
        }
    }
    std::cout << "\nArrays qIt and introIt are equal: " << equal << "\n";
    equal = true;
    for (int i = 0; i < size; i++) {
        if (arrcpy1[i] != arrcpy3[i]) {
            equal = false;
        }
    }
    std::cout << "\nArrays qIt and introRec are equal: " << equal << "\n";

    std::cout << "\n\n";
    */
}





//Block Merge Sort
void mergeBlock(int* arr, int* bfarr, int start, int index, int end) {
    int leftArrayIterator = start;
    int rightArrayIterator = index;
    int bfIterator = start;
    while (leftArrayIterator < index && rightArrayIterator <= end) {
        if (arr[leftArrayIterator] < arr[rightArrayIterator]) {
            bfarr[bfIterator] = arr[leftArrayIterator];
            ++leftArrayIterator;
            ++bfIterator;
        }
        else {
            bfarr[bfIterator] = arr[rightArrayIterator];
            ++rightArrayIterator;
            ++bfIterator;
        }
    }
    for (int k = leftArrayIterator; k < index; k++) {
        bfarr[bfIterator] = arr[k];
        ++bfIterator;
    }
    for (int k = rightArrayIterator; k <= end; k++) {
        bfarr[bfIterator] = arr[k];
        ++bfIterator;
    }
    for (int k = start; k <= end; k++) {
        arr[k] = bfarr[k];
    }
}

void compareHeapRecursiveToIterative(int* arr, int size) {
    int* arrcpy1 = new int[size];
    int* arrcpy2 = new int[size];

    for (int i = 0; i < size; i++) {
        arrcpy1[i] = arr[i];
        arrcpy2[i] = arr[i];
    }


    std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();

    heapSortIntro(arrcpy1, 0, size - 1);

    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
    std::cout << "Recursive HeapSort:\n Elapsed nanoseconds:"
        << std::chrono::duration_cast<std::chrono::nanoseconds> (t1 - t0).count()
        << std::endl;



    std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();

    heapSortIntroIterative(arrcpy2, 0, size - 1);

    std::chrono::steady_clock::time_point t3 = std::chrono::steady_clock::now();
    std::cout << "Iterative HeapSort:\n Elapsed nanoseconds:"
        << std::chrono::duration_cast<std::chrono::nanoseconds> (t3 - t2).count()
        << std::endl;
}



void compareIntroRecursiveToIterative(int* arr, int size) {
    int* arrcpy1 = new int[size];
    int* arrcpy2 = new int[size];

    for (int i = 0; i < size; i++) {
        arrcpy1[i] = arr[i];
        arrcpy2[i] = arr[i];
    }


    std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();

    introSort(arrcpy1, size);

    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
    std::cout << "Recursive IntroSort:\n Elapsed nanoseconds:"
        << std::chrono::duration_cast<std::chrono::nanoseconds> (t1 - t0).count()
        << std::endl;



    std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();

    introSortIterativeQueueLong(arr, size);

    std::chrono::steady_clock::time_point t3 = std::chrono::steady_clock::now();
    std::cout << "Iterative IntroSort:\n Elapsed nanoseconds:"
        << std::chrono::duration_cast<std::chrono::nanoseconds> (t3 - t2).count()
        << std::endl;
}

void bubbleSort(int* arr, int n) {
    bool inversion = true;
    while (inversion) {
        inversion = false;
        for (int i = 1; i < n; i++) {
            if (arr[i - 1] > arr[i]) {
                arr[i] ^= arr[i - 1];
                arr[i - 1] ^= arr[i];
                arr[i] ^= arr[i - 1];
                inversion = true;
            }
        }
    }
}

void selectionSort(int* arr, int n) {
    int tmp;
    for (int i = 0; i < n - 1; i++) {
        for (int j = i + 1; j < n; j++) {
            if (arr[i] > arr[j]) {
                tmp = arr[i];
                arr[i] = arr[j];
                arr[j] = tmp;
            }
        }
    }
}

void insertionSort(int* arr, int n) {
    int value;
    int j;
    for (int i = 1; i < n; i++) {
        value = arr[i];
        j = i;
        while (value < arr[j - 1] && j>0) {
            arr[j] = arr[j - 1];
            --j;
        }
        arr[j] = value;
    }
}

void combSort(int* arr, int n) {
    int gap = n;
    float shrink = 1.3;
    bool swap = true;
    while (swap) {
        gap = gap / shrink;
        if (gap <= 1) {
            gap = 1;
            swap = false;
        }
        int i = 0;
        while (i + gap < n) {
            if (arr[i] > arr[i + gap]) {
                std::swap(arr[i], arr[i + gap]);
                swap = true;
            }
            i ++;
        }
    }
}

void cocktailSort(int* arr, int n) {
    int start = 1;
    int end = n - 1;
    bool swap = true;
    while (swap) {
        int lastswapBegin = start;
        int lastswapEnd = end;
        swap = false;
        for (int i = start; i <= end; i++) {
            if (arr[i] < arr[i - 1]) {
                std::swap(arr[i], arr[i - 1]);
                lastswapEnd = i;
            }
        }
        end = lastswapEnd - 1;
        for (int i = end; i >= start; i--) {
            if (arr[i] < arr[i - 1]) {
                std::swap(arr[i], arr[i - 1]);
                lastswapBegin = i;
                swap = true;
            }
        }
        start = lastswapBegin + 1;
    }
}

void cocktailSortModified(int* arr, int n) {
    int start = 0;
    int end = n - 1;
    int min, minIndex;
    int max, maxIndex;
    while (start<end) {
        min = arr[start];
        minIndex = start;
        max = arr[end];
        maxIndex = end;
        for (int i = start; i <= end; ++i) {
            if (arr[i] > max) {
                max = arr[i];
                maxIndex = i;
            }
            if (arr[i] < min) {
                min = arr[i];
                minIndex = i;
            }
        }
        if (maxIndex == start) {
            maxIndex = minIndex;
        }
        std::swap(arr[start], arr[minIndex]);
        std::swap(arr[end], arr[maxIndex]);
        arr[start] = min;
        arr[end] = max;
        ++start;
        --end;
    }
}









/*
void blockMergeInPlaceWithRotation(int* arr, int start, int index, int end) {
    int positionA = start, lastPositionA, positionB, startingValue;
    printArr(arr, end + 1);
    //std::cout << "\n\n";
    while (positionA < index - 1) {
        lastPositionA = positionA;
        positionA = blockBinarySearchLastIndex(arr, start, index - 1, arr[positionA]);
        if (positionA == index) {
            --positionA;
        }
        positionB = blockBinarySearchForRotation(arr, index, end, arr[positionA]);
        Rotate(arr, lastPositionA, positionB, lastPositionA - positionA);
        printArr(arr, end + 1);
        std::cout << "\n\n";
    }
}
*/






//////////////////////////////////////////////////////////BLOCK MERGE SORT
void Reverse(int* arr, int start, int end) {
    while (start < end) {
        std::swap(arr[start], arr[end]);
        start++, end--;
    }
}

void Rotate(int* arr, int start, int end, int amount=0) {
    if (amount >= 0) {
        if (amount != end - start) {
            Reverse(arr, start, end);
            Reverse(arr, start, start + amount);
            Reverse(arr, start + amount + 1, end);
        }
        return;
    }
    else{
        if (amount != start - end-1) {
            Reverse(arr, start, end);
            Reverse(arr, start, end + amount);
            Reverse(arr, end + amount + 1, end);
        }
        return;
    }
}

/*
void Rotate(int* arr, int start, int end, int amount) {
    amount =
    if (amount >= 0) {
        Reverse(arr, start, end);
        Reverse(arr, start, start + amount);
        Reverse(arr, start + amount + 1, end);
    }
    if (amount < 0) {
        Reverse(arr, start, end);
        Reverse(arr, start, end + amount);
    }
}
*/

void insertionBlockSort(int* arr, int start, int end) {
    int value;
    int j;
    for (int i = start + 1; i <= end; i++) {
        value = arr[i];
        j = i;
        while (value<arr[j - 1] && j>start) {
            arr[j] = arr[j - 1];
            j--;
        }
        arr[j] = value;
    }
}

bool blockBinarySearchElementExists(int* arr, int start, int end, int element) {
    int index;
    while (start <= end) {
        index = start + (end - start) / 2;                     // "advanced" version to avoid overflow ;)
        if (arr[index] == element) {
            return true;
        }
        else {
            if (arr[index] < element) {
                start = index + 1;
            }
            else {
                end = index - 1;
            }
        }
    }
    return false;
}

int blockBinarySearchFirstIndex(int* arr, int start, int end, int element) {
    int index, reservedEnd=end;
    while (start<=end) {
        index = start+(end-start)/ 2;                     // "advanced" version to avoid overflow ;)
        if (arr[index] < element) {
           start = index + 1;
        }
        else {
            end = index - 1;
        }
    }
    if (start>reservedEnd) {
        return reservedEnd;
    }
    return start;
}

int blockBinarySearchLastIndex(int* arr, int start, int end, int element) {
    int index, reservedStart=start;
    while (start <= end) {
        index = start + (end - start) / 2;                  // "advanced" version to avoid overflow ;)
        if (arr[index] > element) {
            end = index-1;
        }
        else {
            start = index + 1;
        }
    }
    if (end<reservedStart) {
        return reservedStart;
    }
    return end;
}

int blockLinearSearchIndex(int* arr, int start, int end, int element) {
    while (arr[start] == element && start < end) {
        ++start;
    }
    return start;
}

bool blockSwap(int* arr, int A_BlockStart, int A_BlockEnd, int B_BlockStart, int B_BlockEnd) {
    int tmp;
    while (A_BlockStart <= A_BlockEnd&&B_BlockStart<=B_BlockEnd) {
        tmp = arr[A_BlockStart];
        arr[A_BlockStart] = arr[B_BlockStart];
        arr[B_BlockStart] = tmp;
        ++A_BlockStart;
        ++B_BlockStart;
    }
    return A_BlockStart > A_BlockEnd;
}

void Redistribute(int* arr, int start, int end, int bufferA, int bufferB) {
    int index = start + bufferA - 1;
    int tmp;
    for (; index >= 0; index--) {
        int i = index;
        tmp = arr[i];
        while (tmp > arr[i + 1] && i < end - bufferB) {
            arr[i] = arr[i + 1];
            ++i;
        }
        arr[i] = tmp;
    }
    index = end-bufferB+1;
    for (; index <=end; index++) {
        int i = index;
        tmp = arr[i];
        while (tmp < arr[i - 1] && i>start) {
            arr[i] = arr[i - 1];
            --i;
        }
        arr[i] = tmp;
    }
}


void blockMergeInternal(int* arr, int A_BlockStart, int A_BlockEnd, int B_BlockStart, int B_BlockEnd, int bufferStart, int bufferEnd) {
    if (B_BlockEnd - B_BlockStart >= 0) {
        int A_swap = 0, B_swap = 0, insert = 0;
        blockSwap(arr, A_BlockStart, A_BlockEnd, bufferStart, bufferEnd);
        while (A_BlockStart + A_swap <= A_BlockEnd && B_BlockStart + B_swap <= B_BlockEnd) {
            if (arr[bufferStart + A_swap] <= arr[B_BlockStart + B_swap]) {
                std::swap(arr[A_BlockStart + insert], arr[bufferStart + A_swap]);
                ++A_swap;
            }
            else {
                std::swap(arr[B_BlockStart + B_swap], arr[A_BlockStart + insert]);
                ++B_swap;
            }
            ++insert;
        }

        for (; A_BlockStart + A_swap <= A_BlockEnd; ++A_swap, ++insert) {
            std::swap(arr[bufferStart + A_swap], arr[A_BlockStart + insert]);
        }
    }
}


/*
void blockMergeInPlace(int* arr, int A_BlockStart, int A_BlockEnd, int B_BlockStart, int B_BlockEnd) {
    int rotationIndex, positionA, amount;
    while(A_BlockEnd>=A_BlockStart&&B_BlockEnd>=B_BlockStart){
        rotationIndex = blockBinarySearchFirstIndex(arr, B_BlockStart, B_BlockEnd, arr[A_BlockStart]);
        if (arr[A_BlockStart] <= arr[rotationIndex]) {
            --rotationIndex;
        }
        else {
            rotationIndex = blockBinarySearchLastIndex(arr, B_BlockStart, B_BlockEnd, arr[rotationIndex]);
        }

        positionA = blockBinarySearchLastIndex(arr, A_BlockStart, A_BlockEnd, arr[A_BlockStart]);
        amount = A_BlockStart - positionA - 1;
        Rotate(arr, A_BlockStart, rotationIndex, amount);
        B_BlockStart += amount;
        A_BlockEnd = B_BlockStart-1;
    }
}
*/

void blockMergeInPlace(int* arr, int A_BlockStart, int A_BlockEnd, int B_BlockStart, int B_BlockEnd) {
    int position;
    while (A_BlockEnd >= A_BlockStart) {
        position = blockBinarySearchFirstIndex(arr, B_BlockStart, B_BlockEnd, arr[A_BlockEnd]);
        if (arr[A_BlockEnd] <= arr[position]) {
            --position;
        }
        Rotate(arr, A_BlockEnd, position, -1);
        --A_BlockEnd;
        --B_BlockStart;
    }
}


void LinearSearchMinTag(int* arr, int start, int end, int blockSize, int& A_MinBlockStart, int& A_MinBlockEnd) {
    int minTag = arr[start + 1];
    int minIndexTag = start;
    for (int i = start + blockSize; i < end; i += blockSize) {
        if (arr[i + 1]<minTag) {
            minTag = arr[i + 1];
            minIndexTag = i;
        }
    }

    A_MinBlockStart = minIndexTag;
    A_MinBlockEnd = minIndexTag + blockSize - 1;
}

void TagA(int* arr, int& start, int& index, int& blockSize, int& bufferA, int& tags) {
    int bufferAIterator = start;
    int startFirstBlockA = start + bufferA;
    int blockAIterator = (startFirstBlockA)+((index - (startFirstBlockA)) % blockSize);
    if (blockAIterator == startFirstBlockA) {
        blockAIterator += blockSize;
    }
    while (bufferAIterator < start + bufferA && blockAIterator < index - 1) {
        std::swap(arr[bufferAIterator], arr[blockAIterator + 1]);
        ++bufferAIterator;
        ++tags;
        blockAIterator += blockSize;
    }
}

void printArr(int* arr, int start, int end) {
    while (start <= end) {
        std::cout << arr[start] << " ";
        ++start;
    }
}


void RollAndDrop(int* arr, int& start, int& index, int& end, int& blockSize, int& bufferA, int& bufferB, int& tags) {
    int size;
    int distance = blockSize - 1;


    int internalEnd = end - bufferB;
    int internalStart = start + bufferA;

    int A_Start = (internalStart)+((index - internalStart) % blockSize);
    int A_End = index-1;

    int A_LastBlockStart = internalStart;
    int A_LastBlockEnd = A_Start - 1;
    if (A_LastBlockEnd < A_LastBlockStart) {
        A_LastBlockEnd = A_LastBlockStart + blockSize - 1;
        A_Start += blockSize;
    }

    int B_Start = index;
    int B_End = end - bufferB;

    int B_LastBlockStart = 0;
    int B_LastBlockEnd = -1;
    int RotationIndex;

    int bufferIterator = start;

    int minA_Start = A_Start;
    int minA_End = minA_Start + distance;


    while (true) {
        if (((B_LastBlockStart <= B_LastBlockEnd) && (arr[minA_Start] <= arr[B_LastBlockEnd])) ||
            (B_Start > B_End)) {
            RotationIndex = blockBinarySearchFirstIndex(arr, B_LastBlockStart, B_LastBlockEnd, arr[minA_Start]);
            if (arr[minA_Start] > arr[RotationIndex]) {
                ++RotationIndex;
            }
            blockSwap(arr, A_Start, A_End, minA_Start, minA_End);
            if (tags>0) {
                std::swap(arr[A_Start + 1], arr[bufferIterator]);
                ++bufferIterator;
                tags--;
            }

            if (RotationIndex != A_Start) {
                Rotate(arr, RotationIndex, A_Start + distance, RotationIndex-A_Start);
            }

            size = B_LastBlockEnd - RotationIndex + 1;
            if (bufferB > 0) {
                blockMergeInternal(arr, A_LastBlockStart, A_LastBlockEnd, A_LastBlockEnd+1, RotationIndex-1, internalEnd + 1, end);
            }
            else {
                blockMergeInPlace(arr, A_LastBlockStart, A_LastBlockEnd, A_LastBlockEnd+1, RotationIndex-1);
            }


            A_LastBlockStart = A_Start - size;
            A_LastBlockEnd = A_LastBlockStart + distance;
            B_LastBlockStart = A_LastBlockEnd + 1;
            B_LastBlockEnd = B_LastBlockStart + size-1;


            A_Start += blockSize;
            if (A_Start > A_End) {
                break;
            }

            LinearSearchMinTag(arr, A_Start, A_End, blockSize, minA_Start, minA_End);
        }
        else if (B_End - B_Start < blockSize) {
            Rotate(arr, A_Start, B_End, B_End - B_Start);
            size = B_End - B_Start + 1;

            B_LastBlockStart = A_Start;
            B_LastBlockEnd = B_LastBlockStart + size-1;
            A_Start += size;
            A_End += size;
            minA_Start += size;
            minA_End += size;
            B_Start += blockSize;
        }
        else {
            blockSwap(arr, A_Start, A_Start + distance, B_Start, B_End);
            B_LastBlockStart = A_Start;
            B_LastBlockEnd = A_Start + distance;
            if (minA_Start == A_Start) {
                minA_Start = A_End + 1;
                minA_End = minA_Start + distance;
            }

            A_Start += blockSize;
            A_End += blockSize;
            B_Start += blockSize;
        }
    }

    if (bufferB > 0) {
        blockMergeInternal(arr, A_LastBlockStart, A_LastBlockEnd, A_LastBlockEnd + 1, internalEnd, internalEnd + 1, end);
    }
    else {
        blockMergeInPlace(arr, A_LastBlockStart, A_LastBlockEnd, A_LastBlockEnd + 1, internalEnd);
    }
}

void blockMerge(int* arr, int start, int index, int end) {
    int totalBufferSize = 2*(int)sqrt(index - start); 
    int bufferA = 0, bufferB=0;
    int startingValue;
    int elementPosition = start;
    int blockSize;
    int integerStep = index - start;


    //CREATE INTERNAL BUFFER (A)
    while (elementPosition < index && bufferA<totalBufferSize/2) {
        startingValue = arr[elementPosition];
        Rotate(arr, start+bufferA, elementPosition);
        ++bufferA;
        elementPosition = blockBinarySearchLastIndex(arr, start+bufferA, index - 1, startingValue);
        if (arr[elementPosition]<=startingValue){
            ++elementPosition;
        }
    }

    
    if (bufferA < 2) {
        blockMergeInPlace(arr, start, index - 1, index, end);
        return;
    }
    

    if(bufferA==totalBufferSize/2){
        //RESET THE POSITION FOR BUFFER (B)
        elementPosition = end;

        //CREATE INTERNAL BUFFER (B)
        while (elementPosition >= index && bufferB < totalBufferSize/2) {
            startingValue = arr[elementPosition];
            if (!blockBinarySearchElementExists(arr,start,start+bufferA-1, startingValue)) {
                Rotate(arr, elementPosition, end-bufferB, -1);
                ++bufferB;
            }
            elementPosition = blockBinarySearchFirstIndex(arr, index, end-bufferB, startingValue);
            if (arr[elementPosition]>=startingValue) {
                --elementPosition;
            }
        }


        //ROTATE INTERNAL BUFFER (B) BACK TO PLACE WHEN THE (B) SIZE IS NOT ENOUGH
        if (bufferA + bufferB < totalBufferSize) {
            while (bufferB > 0) {
                elementPosition = blockBinarySearchFirstIndex(arr, index, end - bufferB, arr[end - bufferB + 1]);
                if (arr[elementPosition] <= arr[end - bufferB + 1]) {
                    elementPosition = blockBinarySearchLastIndex(arr, index, end - bufferB, arr[end - bufferB + 1]);
                    ++elementPosition;
                }
                Rotate(arr, elementPosition, end - bufferB + 1);
                --bufferB;
            }
        }
    }

    if (bufferA + bufferB == totalBufferSize) {
        blockSize = bufferA;
        int tags = 0;
        TagA(arr, start, index, blockSize, bufferA, tags);
        RollAndDrop(arr, start, index, end, blockSize, bufferA, bufferB, tags);
        Redistribute(arr, start, end, bufferA, bufferB);
    }
    else {
        int tags = 0;
        blockSize = integerStep / bufferA + 1;
        TagA(arr, start, index, blockSize, bufferA, tags);
        RollAndDrop(arr, start, index, end, blockSize, bufferA, bufferB, tags);
        Redistribute(arr, start, end, bufferA, bufferB);
    }
}

int blockBinarySearchForRotation(int* arr, int start, int end, int element) {
    int index;
    while (start < end) {
        index = start + (end - start) / 2;
        if (arr[index] > element) {
            start = start + 1;
        }
        else {
            end = end - 1;
        }
    }
    --start;
    return start;
}

int powerOfTwo(int size) {
    int logNum = log2(size);
    return pow(2, logNum);
}

void blockMergeSort(int* arr, int size) {
    float start, index, end;
    int power = powerOfTwo(size);
    float scale = (float)size / (float)power;
    std::cout << scale << "\n";


    for (int merge = 0; merge < power; merge += 16) {
        std::cout << merge << "\n";
        start = (int)(merge * scale);
        end = min((int)((merge + 16) * scale - 1), size - 1);
        std::cout << start << ", " << end << "\n";
        insertionBlockSort(arr, (int)start, (int)end);
    }

    for (int length = 16; length < power; length += length) {
        for (int merge = 0; merge < power; merge += 2 * length) {
            start = (int)(merge * scale);
            index = (int)((merge + length) * scale);
            end = min((int)((merge + length * 2) * scale - 1), size - 1);
            if (arr[(int)end] < arr[(int)start]) {
                Rotate(arr, (int)start, (int)end, (int)(end - index));
            }
            else if (arr[(int)(index - 1)] > arr[(int)index]) {
                blockMerge(arr, (int)start, (int)index, (int)end);
            }
        }
    }
}


void fun(int* arr, int n) {
    int index = n / 2;
    insertionBlockSort(arr, 0, index - 1);
    insertionBlockSort(arr, index, n - 1);
    blockMerge(arr, 0, index, n - 1);
}

void funRand() {
    int n=32;
    int* arr2 = new int[n];
    for (int i = 0; i < n; i++) arr2[i] = rand() % 25;
    int distance = 16;
    for (int i = 0; i < n; i += distance) {
        insertionBlockSort(arr2, i, min(i + distance-1, n - 1));
    }
    printArr(arr2, n);
    blockMerge(arr2, 0, distance, 2*distance-1); 
    printArr(arr2, n);
}

void compareSTLtoSort(int* arr, int n) {
    int* arrcpy1 = new int[n];
    std::vector<int> arrcpy2(n);

    for (int i = 0; i < n; i++) {
        arrcpy1[i] = arr[i];
        arrcpy2[i] = arr[i];
    }

    std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();
    std::sort(arrcpy2.begin(), arrcpy2.end());

    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
    std::cout << "\nstl:\n Elapsed nanoseconds:"
        << std::chrono::duration_cast<std::chrono::nanoseconds> (t1 - t0).count()
        << std::endl;


    std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
    blockMergeSort(arrcpy1, n);
    //fun(arrcpy1, n);


    std::chrono::steady_clock::time_point t3 = std::chrono::steady_clock::now();
    std::cout << "\nblock merge:\n Elapsed nanoseconds:"
        << std::chrono::duration_cast<std::chrono::nanoseconds> (t3 - t2).count()
        << std::endl;

    bool equal = true;
    for (int i = 0; i < n; i++) {
        //std::cout << arrcpy1[i] << " ";
        if (arrcpy1[i] != arrcpy2[i]) {
            //printArr(arrcpy1, n);
            //for (int j = 0; j < n; j++) {
            //    std::cout << arrcpy2[j] << " ";
            //}
            printArr(arr, n);
            std::cout << "\n\nBLOCK MERGE SORT:\n";
            printArr(arrcpy1, n);
            std::cout << "\n\nSTL SORT:\n";
            for (int k = 0; k < n; k++) {
                std::cout << arrcpy2[k] << " ";
            }
            std::cout << "\n\n";
            std::cout << "\narrcpy1 = " << arrcpy1[i] << "\narrcpy2 = " << arrcpy2[i] << " at "<<i<<"\n\n";
            equal = false;
            break;
        }
    }
    assert(equal);
}

int main() {
    int k =1;
    srand(time(NULL));
    for (int j = 0; j < k; j++) {
        int n = 55;//rand()%1000+10000000;//int n = rand()%200;
        std::cout << "# Elements = " << n << "\n";
        //int n = rand()%1000;
        //int arr[] = { 2,1,1,4,2,6,7,0,1,6,6,2,8,1,3,8,0,6,8,8,8,5,8,5,8,8,24,24,24,24,24,22,122,302,102,201,110,100 };
        //int* arr = new int[n]; for (int i = 0; i < n; i++) arr[i] = rand() % 2500;

        int* arr = new int[n]; for (int i = 0; i < n; i++) { 
            arr[i] = rand() % 25000; 
        }

        //int arr[] = { 0, 6, 7, 1, 1, 2, 9, 5, 2, 4, 7, 2, 5, 9, 4, 6, 6, 0, 4, 3, 2, 6, 7, 8, 9, 2, 6, 6, 0, 9, 2, 5, 5, 9, 7, 9, 6, 3 };
        //int arr[] = { 6, 7, 2, 0, 7, 7, 3, 8, 6, 8, 9, 8, 1, 3, 9, 4, 7, 9, 2, 8, 4, 5, 3, 1, 7, 0, 8, 3, 1, 2, 9, 1, 9, 6, 0, 2, 9, 9};
        //int arr[] = { 2, 7, 8, 3, 2, 9, 5, 5, 5, 1, 2, 7, 5, 9, 5, 7, 6, 6, 9, 4, 6, 1, 8, 8, 8, 1, 3, 5, 2, 9, 0, 0, 5, 2, 9, 0, 8, 3, 3, 7, 3, 0, 4, 0, 5, 1, 9, 3, 6, 2, 0, 1 };
        //int arr[] = { 3,0,3,4,4,4,3,1,2,0,0,1,0,0,1,2,1,2,0,1,0,2,3,2,0,2,1,4,0,2,2,1,1,2,3,0,4,4,3,4,3,0,0,1,2,3,4,3,3,2,0,0,1,1,2,3,3,1,0,0,2,2,4,1,1,4,2,3,3,2,3,2,3,2,2,0,1,4,2,3,4,4,1,4,0,0,0 };
        //int arr[] = { 0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,4,4 };

        //FIX THIS!!!!
        
        /*
        int arr[] = { 16,4,1,19,16,20,4,4,13,12,4,7,21,8,2,13,7,15,8,12,9,14,16,14,11,16,9,5,10,16,22,18,12,15,13,1,
            13,12,0,4,3,1,8,4,24,2,6,15,2,15,10,17,4,23,17,2,8,14,23,22,11,15,7,5,13,16,7,18,20,4,1,18,22,19,20,8,24,
            18,24,13,5,8,10,24,11,24,21,3,8,5,22,2,22,20,7,23,5,11,2,17,12,18,21,21,8,19,7,20,11,14,2,18,3,23,19,24,19,
            15,1,10,2,16,23,4,13,3,24,11,20,15,23,7,6,15,24,6,14,15,10,14,21,21,15,5,13,10,7,18,15,1,18,9,8,6,10,2,18,10,
            2,0,5,14,18,16,0,16,14,8,1,7,3,1,11,10,3,24,13,9,18,9,0,1,4,20,6,6,12,15,10,7,14,2,8,2,1,17,9,1,17,22,13,14,19,
            1,8,1,24,15,14,7,22,8,11,3,2,23,15,19,2,6,14,16,5,7,11,9,13,2,11,9,13,24,9,20,0,3,24,5,0,7,10,7,17,9,12,22,20,13,19,6,19,22 };
        */

        //int n = 42;
        //int arr[] = { 11, 9, 4, 5, 13, 0, 19, 17, 11, 1, 1, 7, 17, 16, 16, 9, 15, 17, 18, 12, 0, 13, 10, 16, 9, 17, 7, 5, 7, 14, 12, 2, 5, 6, 4, 19, 17, 16, 3, 16, 13, 10, 0, 8, 13, 19, 12, 16, 7, 15, 10, 17, 18, 15, 13, 15, 7, 4, 0, 2, 16, 8, 0, 4, 6, 18, 7, 9, 14, 5, 14, 8, 7, 9, 8, 11, 18, 9, 18, 8, 10, 4, 6 }; 
        //int arr[] = { 10,0,8,13,19,12,16,7,15,10,17,18,15,13,15,7,4,0,2,16,8,0,4,6,18,7,9,14,5,14,8,7,9,8,11,18,9,18,8,10,4,6 };
        
        //int arr[] = { 15,17,20,21,    2,2,3,4,4,4,5,5,5,5,5,5,5,6,6,7,8,9,10,10,11,   4,5,6,7 };
        //blockMergeInternal(arr, 0, 3, 4, 24, 25, 28);
        //printArr(arr, 29);
        compareSTLtoSort(arr, n);
        //int arr2[] = { 0, 8, 5, 0, 7, 8, 6, 3, 8, 5, 3, 3, 1, 9, 4, 5, 0, 2, 6, 7, 5, 8, 2, 9, 2, 9, 0, 6, 8, 7, 0, 9, 3, 5, 9 };
    }
}