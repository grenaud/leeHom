#ifndef leeHom_Pipeline_h
#define leeHom_Pipeline_h

#include <map>
#include <deque>
#include <mutex>
#include <condition_variable>
#include <stdint.h>

/*!
  Hands chunks of data from a reader to several workers and the results to a writer, in the order
  in which the chunks were submitted whatever the order in which the workers finish them.

     reader --submit()--> [work queue] --take()--> workers --put()--> [reorder buffer] --nextOutput()--> writer

  The chunks are numbered by submit(). nextOutput() returns the result of chunk 0, then 1, then 2 and so on,
  so what the writer produces is in the same order as what the reader read, independently of the number
  of workers.

  At most maxInFlight chunks can have been submitted and not yet been written: submit() waits
  until the writer has caught up, this bounds the memory used.

  Everything blocks on condition variables, nothing polls.
*/
template<class In, class Out>
class OrderedPipeline {
private:
    std::mutex              m;
    std::condition_variable cvWork;    //work available or input closed
    std::condition_variable cvOut;     //a result is available
    std::condition_variable cvSpace;   //the writer made some room

    std::deque< std::pair<uint64_t,In*> > work;
    std::map<uint64_t,Out*>               done;
    uint64_t nextSubmit;
    uint64_t nextWrite;
    uint64_t maxInFlight;
    bool     closed;

public:
    explicit OrderedPipeline(uint64_t maxInFlight_):
	nextSubmit(0), nextWrite(0), maxInFlight(maxInFlight_ < 1 ? 1 : maxInFlight_), closed(false){
    }

    //! Reader: queues a chunk, waits if too many chunks are in flight. The pipeline takes ownership.
    void submit(In * chunk){
	std::unique_lock<std::mutex> lock(m);
	cvSpace.wait(lock, [this]{ return (nextSubmit - nextWrite) < maxInFlight; });
	work.push_back( std::make_pair(nextSubmit, chunk) );
	nextSubmit++;
	cvWork.notify_one();
    }

    //! Reader: no more chunks will be submitted
    void closeInput(){
	std::unique_lock<std::mutex> lock(m);
	closed = true;
	cvWork.notify_all();
	cvOut.notify_all();
    }

    //! Worker: gets the oldest chunk not yet taken. Returns false when there is nothing left to do.
    bool take(uint64_t & rank, In * & chunk){
	std::unique_lock<std::mutex> lock(m);
	cvWork.wait(lock, [this]{ return !work.empty() || closed; });
	if(work.empty())
	    return false;
	rank  = work.front().first;
	chunk = work.front().second;
	work.pop_front();
	return true;
    }

    //! Worker: hands in the result for the chunk numbered rank
    void put(uint64_t rank, Out * out){
	std::unique_lock<std::mutex> lock(m);
	done[rank] = out;
	if(rank == nextWrite)
	    cvOut.notify_one();
    }

    //! Writer: the result of the next chunk in order. Returns false once everything has been written.
    bool nextOutput(Out * & out){
	std::unique_lock<std::mutex> lock(m);
	cvOut.wait(lock, [this]{ return done.count(nextWrite) != 0 || (closed && nextWrite == nextSubmit); });
	typename std::map<uint64_t,Out*>::iterator it = done.find(nextWrite);
	if(it == done.end())
	    return false;
	out = it->second;
	done.erase(it);
	return true;
    }

    //! Writer: the result returned by nextOutput() has been written
    void written(){
	std::unique_lock<std::mutex> lock(m);
	nextWrite++;
	cvSpace.notify_one();
    }
};

#endif
