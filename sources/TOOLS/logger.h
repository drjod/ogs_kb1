#ifndef LOGGER_H
#define LOGGER_H

#include <chrono>
#include <fstream>
#include <cstdio>
#include <string>
#include <sstream>


class Logger
{
	std::string file_name;
	int verbosity;
	bool flag_delete_file;

	Logger() : file_name(""), verbosity(0), flag_delete_file(true)
	{}

	template<typename T>
	std::string collect(T arg)
	{
		std::stringstream ss;
		ss << arg;
		return ss.str();
	}

	template<typename T, typename... Targs>
	std::string collect(T arg, Targs... args)
	{
		std::stringstream ss;
		ss << arg << " " << collect(args...);
		return ss.str();
	}

public:
	static Logger& get_instance()
	{
			static Logger logger;
			return logger;
	}

	//Logger(Logger &logger);// = delete;
	//void operator=(Logger &logger);// = delete;

	void block_deletion()
	{
		flag_delete_file = false;
	}

	void initialize(std::string _dir)
	{
		file_name = _dir + std::string("logger.txt");
	}

	void delete_file()
	{
		if(flag_delete_file)
			std::remove(file_name.c_str());
	}

	void set_verbosity(int _verbosity)
	{
		verbosity = _verbosity;
	}

	template<int VERBOSITY=1, typename... Targs>
	void info(Targs... args)
	{
		if(VERBOSITY <= verbosity)
		{
			std::fstream fout(file_name.c_str(), std::fstream::out | std::fstream:: app);
			std::stringstream ss;

			ss << collect(args...);
			fout << std::string(VERBOSITY-1, '\t') << ss.str() << "\n";
		}
	}

	template<typename... Targs>
	void warning(Targs... args)
	{
		info<1>("WARNING", args...);
	}

};

extern Logger logger;


#endif
