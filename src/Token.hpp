#pragma once

enum class token_status {DATA, BALANCE, END};

class Token
{
public:

  Token()   = default;
  ~Token()  = default;

  virtual bool valid() const final
  {
    return (status_ != token_status::END);
  }

  virtual bool rebalance() const final
  {
    return (status_ == token_status::BALANCE);
  }

  virtual void status(const token_status& s) final
  {
    status_ = s;
  }

  virtual token_status status() const final
  {
    return status_;
  }

  virtual void is_last(const bool b) final
  {
    if (b) {
      status_ = token_status::END;
    }
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
