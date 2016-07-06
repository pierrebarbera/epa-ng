#pragma once

#include "Tiny_Tree.hpp"
#include "Sample.hpp"
#include "MSA_Stream.hpp"

template<class in_type, class out_type>
class Pipeline_Stage {
private:
  Communicator<in_type> acceptor;//acceptor.get(); returns data for next chunk
  Communicator<out_type> forwarder; //forwarder.put(results); sends the results forward
  // MPI that means for example mpi async send to specified address
  // normal that means write output to memory location, or set a pointer
  // possibly collapse into a communication strategy?
  // get and put need to be able to work with different strategies!
  // could work with template?
public:
  void process()// do one chunk, accept, process, forward
  {
    in_type data = acceptor.get();
    out_type results;

    process_(data, results);

    forwarder.put(results);
  }

protected:

  virtual void process_(in_type const& in, out_type& out) = 0;

};

class Placement_Stage : public Pipeline_Stage<MSA_Stream, Sample> {
private:
  std::vector<Tiny_Tree> insertion_trees;
public:
  Placement_Stage (Communicator<MSA_Stream> acceptor, Communicator<Sample> forwarder, std::vector<Tiny_Tree> it)
    : acceptor(acceptor), forwarder(forwarder), insertion_trees(it);

  ~Placement_Stage () = default;

protected:

  void process_(MSA_Stream const& queries, Sample& sample) override
  {
    for (unsigned int local_branch_id = 0; local_branch_id < insertion_trees.size(); ++local_branch_id)
    {
      for (unsigned int cur_seq_id = 0; cur_seq_id < queries.size(); cur_seq_id++)
      {
        sample[cur_seq_id][local_branch_id] = insertion_trees[local_branch_id].place(queries[cur_seq_id]);
      }
    }
  }

};

//Placement_Stage stage(Memory_Communicator<MSA_Stream>(&queries), Memory_Communicator<Sample>(&sample));
