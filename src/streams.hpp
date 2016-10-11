#pragma once

#include <iostream>
#include <string>
#include <streambuf>
#include <memory>

class teebuf : public std::streambuf
{
public:
  teebuf(std::streambuf * sb1, std::streambuf * sb2);
private:
  virtual int overflow(int c);
  virtual int sync();

  std::streambuf * sb1;
  std::streambuf * sb2;
};

class prefixbuf : public std::streambuf
{
public:
  prefixbuf(std::string const& prefix, std::streambuf* sbuf);

private:
  virtual int sync();
  virtual int overflow(int c);

  std::string     prefix;
  std::streambuf* sbuf;
  bool            need_prefix;
};

class Teed_ostream : public std::ostream 
{
public:
  Teed_ostream(std::streambuf * sb1, std::streambuf * sb2);
private:
  teebuf tbuf;
};

class Prefixed_ostream : private virtual prefixbuf, public std::ostream
{
public:
  Prefixed_ostream(std::ostream& out, std::string const& prefix);
};
