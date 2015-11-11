#include "Sequence.h"

#include <string>

using namespace std;

template<typename T, typename D>
Sequence::Sequence(T header, D sequence)
        : sequence(sequence), header(header)
{

}

Sequence::~Sequence()
{

}
