/*
 * log.cpp
 *
 *  Created on: Nov 11, 2009
 *      Author: Adam Auton
 *      ($Revision: 66 $)
 */

#include "output_log.h"

void printLOG(string s)
{
	LOG << s; LOG.flush();
	cout << s; cout.flush();
}

void error(string err_msg, int error_code)
{
	printLOG("Error:" + err_msg + "\n");
	exit(error_code);
}


void error(string err_msg, double value1, double value2, int error_code)
{
	printLOG("Error:" + err_msg + "\n");
	stringstream ss;
	ss << "Value1=" << value1 << " Value2=" << value2 << endl;
	printLOG(ss.str());
	exit(error_code);
}

void warning(string err_msg)
{
	printLOG(err_msg + "\n");
}

void one_off_warning(string err_msg)
{
	static set<string> previous_warnings;
	if (previous_warnings.find(err_msg) == previous_warnings.end())
	{
		printLOG(err_msg + "\n");
		previous_warnings.insert(err_msg);
	}
}

string int2str(int n)
{
  std::ostringstream s2( std::stringstream::out );
  s2 << n;
  return s2.str();
}

string longint2str(long int n)
{
  std::ostringstream s2( std::stringstream::out );
  s2 << n;
  return s2.str();
}

string dbl2str(double n, int prc)
{
  std::ostringstream s2;
  if ( prc > 0 )
    s2.precision(prc);
  s2 << n;
  return s2.str();
}

string dbl2str_fixed(double n, int prc)
{
  std::ostringstream s2;
  s2 << setiosflags( ios::fixed );
  if ( prc > 0 )
    s2.precision(prc);
  s2 << n;
  return s2.str();
}
