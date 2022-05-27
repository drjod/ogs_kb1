#ifndef TIMER_H
#define TIMER_H

#include <chrono>
#include <iostream>


typedef std::chrono::high_resolution_clock _clock;
typedef std::chrono::microseconds precision;

static double total_duration = 0.;

template<typename Stream=std::ostream>
struct Timer
{
	Timer(const char* _identifier="", Stream& _stream=std::cout) : stream(_stream), identifier(_identifier)
	{
		start = _clock::now();
		stream << "Timer " << identifier;
	}
	~Timer()
	{
		const double duration = std::chrono::duration_cast<precision>(_clock::now() - start).count();
		total_duration += duration;
		stream << ":\t" << duration << "\t" << total_duration << '\n';
	}
private:
	Stream& stream;
	const char* identifier;
	std::chrono::time_point<_clock> start;

};

#endif
