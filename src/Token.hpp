#pragma once

class Token
{
public:
  enum class token_status {DATA, BALANCE, END};

  Token()   = default;
  ~Token()  = default;

  virtual operator bool() const final
  {
    return (status_ == token_status::DATA);
  }

  virtual bool rebalance() const final
  {
    return (status_ == token_status::BALANCE);
  }

  virtual void status(token_status s) final
  {
    status_ = s;
  }

private:

  token_status status_ = token_status::DATA;
  
};

class VoidToken : public Token
{
public:
  VoidToken() = default;
  ~VoidToken()= default;
};
