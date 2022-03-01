#ifndef HISTORY_QUEUE_H
#define HISTORY_QUEUE_H

#include <iostream>
#include <vector>

template<class T, size_t MAX_SIZE=20>
class HistoryQueue
{
private:
    std::vector<T> past;
    std::vector<T> future;

public:
    void Clear() { past.clear(); future.clear(); }

    void Insert(T ele)
    {
        // insert element
        past.push_back(ele);
        // clear redo records
        future.clear();
    }

    void Undo()
    {
        if (past.empty())
            return;
        // insert to redo record
        future.push_back(past.back());
        // delete from undo record
        past.pop_back();
    }

    void Redo()
    {
        if (future.empty())
            return;
        // insert to undo record
        past.push_back(future.back());
        // delete from redo record
        future.pop_back();
    }

    const std::vector<T>& HistoryList() const
    {
        return past;
    }
};









#endif