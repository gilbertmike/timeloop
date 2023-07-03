#include <fstream>

#include "util/accelergy_interface.hpp"
#include "util/banner.hpp"

#include "applications/looptree-model/model.hpp"
#include "loop-analysis/nest-analysis.hpp"
#include "loop-analysis/isl-analysis/isl-nest-analysis.hpp"
#include "loop-analysis/mapping-to-isl/fused-mapping-to-isl.hpp"
#include "loop-analysis/isl-ir.hpp"
#include "isl-wrapper/ctx-manager.hpp"
#include "mapping/fused-mapping.hpp"
#include "workload/fused-workload.hpp"
#include <isl/constraint.h>
#include <barvinok/isl.h>
#include "loop-analysis/temporal-analysis.hpp"
#include "loop-analysis/spatial-analysis.hpp"

//--------------------------------------------//
//                Application                 //
//--------------------------------------------//

char* isl_pw_qpolynomial_fold_to_str(isl_pw_qpolynomial_fold* pwqf)
{
  auto p_printer = isl_printer_to_str(GetIslCtx().get());
  p_printer = isl_printer_print_pw_qpolynomial_fold(p_printer, pwqf);
  auto p_str = isl_printer_get_str(p_printer);
  isl_printer_free(p_printer);
  return p_str;
}

char* isl_pw_qpolynomial_fold_to_str(isl_pw_qpolynomial_fold* pwqf, void*)
{
  auto p_printer = isl_printer_to_str(GetIslCtx().get());
  p_printer = isl_printer_print_pw_qpolynomial_fold(p_printer, pwqf);
  auto p_str = isl_printer_get_str(p_printer);
  isl_printer_free(p_printer);
  return p_str;
}

isl_stat foo2(isl_qpolynomial* qp, void* qp_out)
{
  *((isl_pw_qpolynomial**) qp_out) = isl_pw_qpolynomial_from_qpolynomial(qp);
  std::cout << isl_pw_qpolynomial_to_str(isl_pw_qpolynomial_from_qpolynomial(
    qp
  )) << std::endl;
  return isl_stat_ok;
}

isl_bool foo(isl_set* set, isl_qpolynomial_fold* fold, void* qp_out)
{
  std::cout << isl_set_to_str(set) << std::endl;
  isl_qpolynomial_fold_foreach_qpolynomial(
    fold,
    foo2,
    qp_out
  );
  return isl_bool_true;
}

template <class Archive>
void Application::serialize(Archive& ar, const unsigned int version)
{
  if (version == 0)
  {
    ar& BOOST_SERIALIZATION_NVP(workload_);
  }
}

Application::Application(config::CompoundConfig* config,
                         std::string output_dir,
                         std::string name) :
    name_(name)
{    
  auto rootNode = config->getRoot();

  // Model application configuration.
  auto_bypass_on_failure_ = false;
  std::string semi_qualified_prefix = name;

  if (rootNode.exists("model"))
  {
    auto model = rootNode.lookup("model");
    model.lookupValue("verbose", verbose_);
    model.lookupValue("auto_bypass_on_failure", auto_bypass_on_failure_);
    model.lookupValue("out_prefix", semi_qualified_prefix);
  }

  out_prefix_ = output_dir + "/" + semi_qualified_prefix;

  if (verbose_)
  {
    for (auto& line: banner)
      std::cout << line << std::endl;
    std::cout << std::endl;
  }

  auto workload = problem::ParseFusedWorkload(rootNode.lookup("problem"));

  // Architecture configuration.
  config::CompoundConfigNode arch;
  if (rootNode.exists("arch"))
  {
    arch = rootNode.lookup("arch");
  }
  else if (rootNode.exists("architecture"))
  {
    arch = rootNode.lookup("architecture");
  }
  
  bool is_sparse_topology = rootNode.exists("sparse_optimizations");
  arch_specs_ = model::Engine::ParseSpecs(arch, is_sparse_topology);

  if (rootNode.exists("ERT"))
  {
    auto ert = rootNode.lookup("ERT");
    if (verbose_)
      std::cout << "Found Accelergy ERT (energy reference table), replacing internal energy model." << std::endl;
    arch_specs_.topology.ParseAccelergyERT(ert);
    if (rootNode.exists("ART")){ // Nellie: well, if the users have the version of Accelergy that generates ART
      auto art = rootNode.lookup("ART");
      if (verbose_)
        std::cout << "Found Accelergy ART (area reference table), replacing internal area model." << std::endl;
      arch_specs_.topology.ParseAccelergyART(art);  
    }
  }
  else
  {
#ifdef USE_ACCELERGY
    // Call accelergy ERT with all input files
    if (arch.exists("subtree") || arch.exists("local"))
    {
      accelergy::invokeAccelergy(config->inFiles, semi_qualified_prefix, output_dir);
      std::string ertPath = out_prefix_ + ".ERT.yaml";
      auto ertConfig = new config::CompoundConfig(ertPath.c_str());
      auto ert = ertConfig->getRoot().lookup("ERT");
      if (verbose_)
        std::cout << "Generate Accelergy ERT (energy reference table) to replace internal energy model." << std::endl;
      arch_specs_.topology.ParseAccelergyERT(ert);
        
      std::string artPath = out_prefix_ + ".ART.yaml";
      auto artConfig = new config::CompoundConfig(artPath.c_str());
      auto art = artConfig->getRoot().lookup("ART");
      if (verbose_)
        std::cout << "Generate Accelergy ART (area reference table) to replace internal area model." << std::endl;
      arch_specs_.topology.ParseAccelergyART(art);
    }
#endif
  }

  // Sparse optimizations
  config::CompoundConfigNode sparse_optimizations;
  if (is_sparse_topology)
    sparse_optimizations = rootNode.lookup("sparse_optimizations");
      sparse_optimizations_ = new sparse::SparseOptimizationInfo(sparse::ParseAndConstruct(sparse_optimizations, arch_specs_));
  // characterize workload on whether it has metadata
  workload_.SetDefaultDenseTensorFlag(sparse_optimizations_->compression_info.all_ranks_default_dense);
  
  if (verbose_)
    std::cout << "Sparse optimization configuration complete." << std::endl;

  arch_props_ = new ArchProperties(arch_specs_);
  // Architecture constraints.
  config::CompoundConfigNode arch_constraints;

  if (arch.exists("constraints"))
    arch_constraints = arch.lookup("constraints");
  else if (rootNode.exists("arch_constraints"))
    arch_constraints = rootNode.lookup("arch_constraints");
  else if (rootNode.exists("architecture_constraints"))
    arch_constraints = rootNode.lookup("architecture_constraints");

  constraints_ = new mapping::Constraints(*arch_props_, workload_);
  constraints_->Parse(arch_constraints);

  if (verbose_)
    std::cout << "Architecture configuration complete." << std::endl;

  mapping::FusedMapping mapping =
    mapping::ParseMapping(rootNode.lookup("mapping"), workload);

  auto mapping_analysis_result =
    analysis::OccupanciesFromMapping(mapping, workload);

  for (const auto& [buf, occ] : mapping_analysis_result.lbuf_to_occupancy)
  {
    auto result = analysis::TemporalReuseAnalysis(
      analysis::TemporalReuseAnalysisInput(
        occ,
        analysis::BufTemporalReuseOpts{
          .exploit_temporal_reuse=1
        }
      )
    );

    auto p_fill = result.fill.map.copy();
    auto p_fill_count = isl_pw_qpolynomial_sum(isl_map_card(p_fill));

    auto p_occ = result.effective_occupancy.map.copy();
    isl_bool tight;
    auto p_occ_count = isl_pw_qpolynomial_bound(
      isl_map_card(p_occ),
      isl_fold_max,
      &tight
    );
    assert(tight == isl_bool_true);

    std::cout << "[Occupancy]" << buf << ": "
      << isl_pw_qpolynomial_fold_to_str(p_occ_count) << std::endl;
    std::cout << "[Fill]" << buf << ": "
      << isl_pw_qpolynomial_to_str(p_fill_count) << std::endl;
    isl_pw_qpolynomial_free(p_fill_count);
  }

  auto p_fold = isl_pw_qpolynomial_fold_read_from_str(
    GetIslCtx().get(),
    "{ [x, y] -> max(25 + 5*x + y) }"
  );
  auto p_map = isl_map_read_from_str(
    GetIslCtx().get(),
    "{ [x] -> [x, y] : 0 <= y < 5 }"
  );
  isl_pw_qpolynomial* p_qp;
  isl_pw_qpolynomial_fold_every_piece(
    p_fold,
    foo,
    &p_qp
  );
  p_qp = isl_map_apply_pw_qpolynomial(p_map, p_qp);
  std::cout << isl_pw_qpolynomial_to_str(p_qp) << std::endl;

  for (const auto& [compute, tiling] : mapping_analysis_result.branch_tiling)
  {
    auto p_ops = isl_map_card(tiling.copy());
    std::cout << "[Operations]" << compute << ": "
      << isl_pw_qpolynomial_to_str(p_ops) << std::endl;

    auto assumed_parallelism =
      mapping_analysis_result.compute_to_assumed_parallelism.at(compute);
    auto p_parallelism = isl_val_int_from_si(
      GetIslCtx().get(),
      static_cast<int>(assumed_parallelism)
    );
    auto p_latency = isl_pw_qpolynomial_scale_down_val(
      isl_pw_qpolynomial_copy(p_ops),
      p_parallelism
    );
    mapping_analysis_result.compute_latency_aggregator.SetLatency(
      compute,
      p_latency
    );
    isl_pw_qpolynomial_free(p_ops);
  }

  mapping_analysis_result.compute_latency_aggregator.CalculateLatency();

  // // for (const auto& [buf, fill] : fills)
  // // {
  // //   std::cout << buf << std::endl;
  // //   auto p_fill_count = isl_map_card(fill.map.copy());
  // //   std::cout << isl_pw_qpolynomial_to_str(p_fill_count) << std::endl;
  // //   isl_pw_qpolynomial_free(p_fill_count);
  // // }

  // for (const auto& [buf, occ] : occupancies)
  // {
  //   std::cout << buf << std::endl;
  //   auto p_occ_count = isl_map_card(occ.map.copy());
  //   std::cout << isl_pw_qpolynomial_to_str(p_occ_count) << std::endl;
  //   isl_pw_qpolynomial_free(p_occ_count);
  // }
}

Application::~Application()
{
  if (arch_props_)
    delete arch_props_;

  if (constraints_)
    delete constraints_;

  if (sparse_optimizations_)
    delete sparse_optimizations_;
}

// Run the evaluation.
Application::Stats Application::Run()
{
  // Output file names.
  model::Engine engine;
  engine.Spec(arch_specs_);

  auto level_names = arch_specs_.topology.LevelNames();

  // if (engine.IsEvaluated())
  // {
  //   if (!sparse_optimizations_->no_optimization_applied)
  //   {   
  //     std::cout << "Utilization = " << std::setw(4) << OUT_FLOAT_FORMAT << std::setprecision(2) << engine.Utilization()
  //             << " | pJ/Algorithmic-Compute = " << std::setw(8) << OUT_FLOAT_FORMAT << PRINTFLOAT_PRECISION << engine.Energy() /
  //     engine.GetTopology().AlgorithmicComputes()
  //             << " | pJ/Compute = " << std::setw(8) << OUT_FLOAT_FORMAT << PRINTFLOAT_PRECISION << engine.Energy() /
  //     engine.GetTopology().ActualComputes() << std::endl;
  //   }
  //   else
  //   {
  //     std::cout << "Utilization = " << std::setw(4) << OUT_FLOAT_FORMAT << std::setprecision(2) << engine.Utilization()
  //                << " | pJ/Compute = " << std::setw(8) << OUT_FLOAT_FORMAT << PRINTFLOAT_PRECISION << engine.Energy() /
  //     engine.GetTopology().ActualComputes() << std::endl;
  //   }
  //   std::ofstream map_txt_file(map_txt_file_name);
  //   mapping.PrettyPrint(map_txt_file, arch_specs_.topology.StorageLevelNames(), engine.GetTopology().UtilizedCapacities(), engine.GetTopology().TileSizes());
  //   map_txt_file.close();

  //   std::ofstream stats_file(stats_file_name);
  //   stats_file << engine << std::endl;
  //   stats_file.close();
  // }

  Stats stats;
  // stats.cycles = engine.Cycles();
  // stats.energy = engine.Energy();
  return stats;
}