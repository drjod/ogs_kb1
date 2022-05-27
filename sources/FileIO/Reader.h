/**
 * \file Reader.h
 * 13/02/2012 LB Initial implementation
 */

#ifndef READER_H
#define READER_H

#include <string>

namespace FileIO
{

/// @brief Base class which enables reading an object from string, stringstream
/// or file.
///
/// When subclassing you only need to implement void read(std::ostream& stream).
class Reader
{
public:
	Reader();
	virtual ~Reader() {};

	/// @brief Reads the object from a string.
	void readFromString(std::string str);

	/// @brief Reads the object from the given file.
	void readFromFile(std::string filename);

protected:
	/// @brief Reads an object from the given stream.
	/// This method must be implemented by a subclass.
	virtual void read(std::istream& stream) = 0;
	
	/// @brief The stream to read from.
	std::stringstream _stream;

private:
	/* data */
};

} // namespace FileIO


#endif // READER_H
