#pragma once

#include <iostream>
#include <iomanip>

#include <map>
#include <string>

#include <chrono>
#include <atomic>

#include "fmmtl/common.hpp"
#include "fmmtl/config.hpp"

/** An RAII class
 * Updates a listener with the amount of time the Ticker was alive.
 */
template <typename Listener>
class TickerNotifier {
 public:
  typedef std::chrono::high_resolution_clock clock;
  typedef typename clock::time_point         time_point;
  typedef typename clock::duration           duration;
  typedef typename duration::rep             tick_type;

  TickerNotifier()
      : owner_(nullptr), starttime_(clock::now()) {}
  explicit TickerNotifier(Listener* owner)
      : owner_(owner), starttime_(clock::now()) {}
  // Allow moving
  TickerNotifier(TickerNotifier&& t) :
      owner_(t.owner_), starttime_(t.starttime_) {
    t.owner_ = nullptr;
  }
  // Disable copying
  TickerNotifier(const TickerNotifier&) = delete;
  TickerNotifier& operator=(const TickerNotifier&) = delete;
  // Destructor
  ~TickerNotifier() {
    duration tick_time = elapsed();
    if (owner_)
      owner_->operator+=(tick_time);
  }
  // Get the duration on this Ticker
  duration elapsed() const {
    return clock::now() - starttime_;
  }
  // Get the seconds on this Ticker
  double seconds() const {
    typedef std::chrono::duration<double> units;
    return std::chrono::duration_cast<units>(elapsed()).count();
  }
 private:
  Listener* owner_;
  time_point starttime_;
};
/** A quick class for timing code:
 * Usage:
 *
 * Ticker ticker;
 * // code to time
 * double time = ticker.seconds();
 */
typedef TickerNotifier<std::chrono::duration<double>> Ticker;


/** A quick class for timing segments of code:
 * Usage:
 *
 * Timer timer;
 * { auto time = timer.start();
 *   // code to time
 * }
 * std::cout << timer << std::endl;
 */
class Timer {
 public:
  typedef TickerNotifier<Timer>       ticker_type;
  typedef ticker_type::clock          clock;
  typedef typename clock::time_point  time_point;
  typedef typename clock::duration    duration;
  typedef typename duration::rep      tick_type;

  // Start by returning an RAII ticker
  ticker_type start() {
    return ticker_type(this);
  }
  // Add a duration
  void operator+=(const duration& d) {
    ticks_ += d.count();
  }
  // Reset this timer
  void reset() {
    ticks_ = tick_type(0);
  }
  // Get the duration on this Timer
  duration total() const {
    return duration(ticks_);
  }
  // Get the seconds on this Timer
  double seconds() const {
    typedef std::chrono::duration<double> units;
    return std::chrono::duration_cast<units>(total()).count();
  }
  // Print this Timer
  friend std::ostream& operator<<(std::ostream& s, const Timer& t) {
    return s << t.seconds() << "secs";
  }
 private:
  //std::atomic<tick_type> ticks_;
  tick_type ticks_;
};


/** @class Logger
 * @brief Logging class to keep the hit count and total time of sections of code
 *
 * Usage:
 * Logger logger;
 *
 * { Logger::timer timer = logger.log("identifier");
 *  // code to track
 * }
 *
 * // Print all string-identified events' hit counts and timings
 * std::cout << log << std::endl;
 */
class Logger {
  struct EventData;

 public:
  typedef TickerNotifier<EventData>   ticker_type;
  typedef ticker_type::clock          clock;
  typedef typename clock::time_point  time_point;
  typedef typename clock::duration    duration;
  typedef typename duration::rep      tick_type;

  /** Start a ticker for an event. */
  inline ticker_type log(const std::string& event) {
    auto range = data_.equal_range(event);
    if (range.first == range.second) {
#pragma omp critical
      {
      range.first = data_.insert(range.first,
                                 std::make_pair(event, new EventData()));
      }
    }
    return ticker_type((*range.first).second);
  }

  /** Erase all events in timer */
  inline void clear() {
#pragma omp critical
    {
      for (auto& event : data_)
        delete event.second;
      data_.clear();
    }
  }

  //! Print all events and timing to an ostream
  friend std::ostream& operator<<(std::ostream& s, const Logger& log) {
    for (auto& event : log.data_)
      s << std::setw(20) << std::left << event.first
        << ": " << *(event.second) << std::endl;
    return s;
  }

 private:
  struct EventData {
    void operator+=(const duration& time) {
      total_ += time;
      ++call_;
    }
    double seconds() const {
      return total_.seconds();
    }
    unsigned calls() const {
      return call_;
    }
    friend std::ostream& operator<<(std::ostream& s, const EventData& e) {
      return s << e.calls() << " (calls) * "
               << e.seconds() / e.calls() << " (sec/call) = "
               << e.seconds() << " (sec)";
    }
   private:
    Timer total_;
    //std::atomic<unsigned> call_; // Not required if threads have unique event strings
    unsigned call_;
  };

  // A map of pointers to (total_time, #calls) with string identifiers
  std::map<std::string, EventData*> data_;
};

//! Global static logger rather than a singleton for efficiency/consistency
static Logger fmmtl_logger;


#if defined(FMMTL_LOGGING)
#define FMMTL_LOG(STRING) auto t##__LINE__ = fmmtl_logger.log(std::string(STRING) + " [" + std::to_string(omp_get_thread_num()) + ']')
#define FMMTL_PRINT_LOG(OUT) OUT << fmmtl_logger << std::endl
#else
#define FMMTL_LOG(STRING)
#define FMMTL_PRINT_LOG(OUT)
#endif
