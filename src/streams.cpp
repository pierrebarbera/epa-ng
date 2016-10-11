#include "streams.hpp"

teebuf::teebuf(std::streambuf * sb1, std::streambuf * sb2)
    : sb1(sb1)
    , sb2(sb2)
  { }

int teebuf::overflow(int c)
{
  if (c == EOF)
  {
    return !EOF;
  }
  else
  {
    int r1 = EOF;
    int r2 = EOF;

    if (sb1)
      r1 = sb1->sputc(c);
    if (sb2)    
      r2 = sb2->sputc(c);

    return r1 == EOF || r2 == EOF ? EOF : c;
  }
}
  
// Sync both teed buffers.
int teebuf::sync()
{
  int r1 = EOF;
  int r2 = EOF;

  if (sb1)
    r1 = sb1->pubsync();
  if (sb2)
    r2 = sb2->pubsync();

  return r1 == 0 && r2 == 0 ? 0 : -1;
}

prefixbuf::prefixbuf(std::string const& prefix, std::streambuf* sbuf)
  : prefix(prefix)
  , sbuf(sbuf)
  , need_prefix(true)
  { }

int prefixbuf::sync() 
{
  return this->sbuf->pubsync();
}

int prefixbuf::overflow(int c) 
{
  if (c != std::char_traits<char>::eof()) 
  {
    if (this->need_prefix
      && !this->prefix.empty()
      && this->prefix.size() != 
        (unsigned int) this->sbuf->sputn(&this->prefix[0], this->prefix.size())) 
    {
      return std::char_traits<char>::eof();
    }
    this->need_prefix = c == '\n';
  }
  return this->sbuf->sputc(c);
}

Teed_ostream::Teed_ostream(std::streambuf * sb1, std::streambuf * sb2)
  : std::ostream(&tbuf)
  , tbuf(sb1, sb2)
  { }

Prefixed_ostream::Prefixed_ostream(std::ostream& out, std::string const& prefix)
  : prefixbuf(prefix, out.rdbuf())
  , std::ios(static_cast<std::streambuf*>(this))
  , std::ostream(static_cast<std::streambuf*>(this)) 
  { }
