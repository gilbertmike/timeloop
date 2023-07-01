#include "barvinok/isl.h"
#include "isl/constraint.h"

#include "isl-wrapper/ctx-manager.hpp"
#include "isl-wrapper/isl-functions.hpp"

namespace isl {

size_t dim(const map& map, isl_dim_type dim_type)
{
  return isl_map_dim(map.get(), dim_type);
}

map dim_projector(space space, size_t start, size_t n)
{
  return isl::manage(dim_projector(
    space.release(),
    start,
    n
  ));
}

isl_map* dim_projector(isl_space* space, size_t start, size_t n)
{
  return isl_map_project_out(
    isl_map_identity(isl_space_map_from_set(space)),
    isl_dim_in,
    start,
    n
  );
}

map project_dim(map map, isl_dim_type dim_type, size_t start, size_t n)
{
  return isl::manage(isl_map_project_out(map.release(), dim_type, start, n));
}

map project_dim_in_after(map map, size_t start)
{
  auto n_dim_in = isl_map_dim(map.get(), isl_dim_in);
  return project_dim(map, isl_dim_in, start, n_dim_in - start);
}

isl::map project_last_dim(isl::map map)
{
  auto n_dim_in = isl_map_dim(map.get(), isl_dim_in);
  return project_dim(map, isl_dim_in, n_dim_in-1, 1);
}

isl::map map_from_multi_aff(isl::multi_aff maff)
{
  return isl::manage(isl_map_from_multi_aff(maff.release()));
}
isl::map map_from_multi_aff(isl::pw_multi_aff maff)
{
  return isl::manage(isl_map_from_pw_multi_aff(maff.release()));
}

isl::space
space_alloc(isl::ctx ctx, size_t n_params, size_t n_dim_in, size_t n_dim_out)
{
  return isl::manage(isl_space_alloc(ctx.release(),
                                     n_params,
                                     n_dim_in,
                                     n_dim_out));
}

isl::space
space_set_alloc(isl::ctx ctx, size_t n_params, size_t n_dim_set)
{
  return isl::manage(isl_space_set_alloc(ctx.release(), n_params, n_dim_set));
}

isl::space
space_from_domain_and_range(isl::space domain, isl::space range)
{
  return isl::manage(isl_space_map_from_domain_and_range(
    domain.release(),
    range.release()
  ));
}

isl::aff
set_coefficient_si(isl::aff aff, isl_dim_type dim_type, size_t pos, int val)
{
  return isl::manage(isl_aff_set_coefficient_si(aff.release(),
                                                dim_type,
                                                pos,
                                                val));
}

isl::aff set_constant_si(isl::aff aff, int val)
{
  return isl::manage(isl_aff_set_constant_si(aff.release(), val));
}

isl::aff si_on_domain(isl::space space, int val)
{
  return isl::manage(isl_aff_val_on_domain_space(
    space.release(),
    isl_val_int_from_si(space.ctx().get(), val)
  ));
}

isl::map add_dims(isl::map map, isl_dim_type dim_type, size_t n_dims)
{
  return isl::manage(isl_map_add_dims(map.release(), dim_type, n_dims));
}

isl::map insert_dims(isl::map map,
                     isl_dim_type dim_type, size_t pos, size_t n_dims)
{
  return isl::manage(isl_map_insert_dims(map.release(),
                                         dim_type, pos, n_dims));
}

isl::map move_dims(isl::map map,
                   isl_dim_type dst_dim_type, size_t dst,
                   isl_dim_type src_dim_type, size_t src,
                   size_t n_dims)
{
  return isl::manage(
    isl_map_move_dims(map.release(),
                      dst_dim_type, dst,
                      src_dim_type, src,
                      n_dims)
  );
}

map map_to_shifted(space domain_space, size_t pos, int shift)
{
  return isl::manage(map_to_shifted(domain_space.release(), pos, shift));
}

isl_map*
map_to_shifted(__isl_take isl_space* domain_space, size_t pos, int shift)
{
  auto p_maff = isl_multi_aff_identity_on_domain_space(domain_space);
  p_maff = isl_multi_aff_set_at(
    p_maff,
    pos,
    isl_aff_set_constant_si(isl_multi_aff_get_at(p_maff, pos), shift)
  );
  return isl_map_from_multi_aff(p_maff);
}

map map_to_all_at_dim(space domain_space, size_t pos)
{
  auto p_domain_space = domain_space.release();
  auto p_map = isl_map_identity(isl_space_map_from_set(p_domain_space));
  p_map = isl_map_project_out(p_map, isl_dim_out, pos, 1);
  p_map = isl_map_insert_dims(p_map, isl_dim_out, pos, 1);
  return isl::manage(p_map);
}

isl::map fix_si(isl::map map, isl_dim_type dim_type, size_t pos, int val)
{
  return isl::manage(isl_map_fix_si(map.release(), dim_type, pos, val));
}

map insert_equal_dims(map map, size_t in_pos, size_t out_pos, size_t n)
{
  return isl::manage(insert_equal_dims(map.release(), in_pos, out_pos, n));
}

isl_map*
insert_equal_dims(__isl_take isl_map* p_map, size_t in_pos, size_t out_pos,
                  size_t n)
{
  p_map = isl_map_insert_dims(p_map, isl_dim_in, in_pos, n);
  p_map = isl_map_insert_dims(p_map, isl_dim_out, out_pos, n);

  auto p_ls = isl_local_space_from_space(isl_map_get_space(p_map));
  for (size_t i = 0; i < n; ++i)
  {
    auto c = isl_constraint_alloc_equality(isl_local_space_copy(p_ls));
    c = isl_constraint_set_coefficient_si(c, isl_dim_in, i + in_pos, 1);
    c = isl_constraint_set_coefficient_si(c, isl_dim_out, i + out_pos, -1);
    p_map = isl_map_add_constraint(p_map, c);
  }
  isl_local_space_free(p_ls);

  return p_map;
}

isl::multi_aff
insert_equal_dims(isl::multi_aff maff, size_t in_pos, size_t out_pos, size_t n)
{
  auto p_maff = maff.release();
  p_maff = isl_multi_aff_insert_dims(p_maff, isl_dim_in, in_pos, n);
  p_maff = isl_multi_aff_insert_dims(p_maff, isl_dim_out, out_pos, n);

  for (size_t i = 0; i < n; ++i)
  {
    auto p_aff = isl_multi_aff_get_at(p_maff, out_pos+i);
    p_aff = isl_aff_set_coefficient_si(p_aff, isl_dim_in, in_pos+i, 1);
    p_maff = isl_multi_aff_set_at(p_maff, out_pos+i, p_aff);
  }

  return isl::manage(p_maff);
}

__isl_give isl_multi_aff*
insert_equal_dims(__isl_take isl_multi_aff* p_maff,
                  int in_pos, int out_pos, int n)
{
  p_maff = isl_multi_aff_insert_dims(p_maff, isl_dim_in, in_pos, n);
  p_maff = isl_multi_aff_insert_dims(p_maff, isl_dim_out, out_pos, n);

  for (auto i = 0; i < n; ++i)
  {
    auto p_aff = isl_multi_aff_get_at(p_maff, out_pos+i);
    p_aff = isl_aff_set_coefficient_si(p_aff, isl_dim_in, in_pos+i, 1);
    p_maff = isl_multi_aff_set_at(p_maff, out_pos+i, p_aff);
  }

  return p_maff;
}

isl::map insert_dummy_dim_ins(isl::map map, size_t pos, size_t n)
{
  auto p_map = map.release();
  p_map = isl_map_insert_dims(p_map, isl_dim_in, pos, n);

  auto p_ls = isl_local_space_from_space(isl_map_get_space(p_map));
  for (size_t i = 0; i < n; ++i)
  {
    auto c = isl_constraint_alloc_equality(isl_local_space_copy(p_ls));
    c = isl_constraint_set_coefficient_si(c, isl_dim_in, i + pos, 1);
    c = isl_constraint_set_constant_si(c, 0);
    p_map = isl_map_add_constraint(p_map, c);
  }
  isl_local_space_free(p_ls);

  return isl::manage(p_map);
}

isl_pw_qpolynomial* sum_map_range_card(isl::map map)
{
  auto p_domain = map.domain().release();
  auto p_count = isl_map_card(map.release());
  return isl_set_apply_pw_qpolynomial(p_domain, p_count);
}

isl_pw_qpolynomial* set_card(isl::set set)
{
  return isl_set_card(set.release());
}

double val_to_double(isl_val* val)
{
  auto num = isl_val_get_num_si(val);
  auto den = isl_val_get_den_si(val);
  isl_val_free(val);
  return (double)num / (double)den;
}

isl_val* get_val_from_singular_qpolynomial(isl_pw_qpolynomial* pw_qp)
{
  return isl_pw_qpolynomial_eval(
    pw_qp,
    isl_point_zero(isl_pw_qpolynomial_get_domain_space(pw_qp))
  );
}

isl_val* get_val_from_singular_qpolynomial_fold(isl_pw_qpolynomial_fold* pwf)
{
  return isl_pw_qpolynomial_fold_eval(
    pwf,
    isl_point_zero(isl_pw_qpolynomial_fold_get_domain_space(pwf))
  );
}

isl::map ConstraintDimEquals(isl::map map, size_t n_dims)
{
  auto p_map = map.release();
  auto p_space = isl_map_get_space(p_map);
  auto p_ls = isl_local_space_from_space(p_space);

  isl_constraint* p_c;
  for (size_t i = 0; i < n_dims; ++i)
  {
    p_c = isl_constraint_alloc_equality(isl_local_space_copy(p_ls));
    p_c = isl_constraint_set_coefficient_si(p_c, isl_dim_in, i, 1);
    p_c = isl_constraint_set_coefficient_si(p_c, isl_dim_out, i, -1);
    p_map = isl_map_add_constraint(p_map, p_c);
  }
  isl_local_space_free(p_ls);

  return isl::manage(p_map);
}

isl::map MapToPriorData(size_t n_in_dims, size_t top)
{
  isl_space* p_space;
  isl_map* p_map;
  isl_local_space* p_ls;
  isl_constraint* p_c;

  // Goal: { [i0, ..., i{n_in_dims-1}] -> [i0, ..., i{top}-1, o{top+1}, ..., o{n_in_dims}] }
  p_space = isl_space_alloc(GetIslCtx().get(), 0, n_in_dims, n_in_dims);
  p_map = isl_map_empty(isl_space_copy(p_space));
  p_ls = isl_local_space_from_space(p_space);

  if (top > 0)
  {
    auto p_tmp_map = isl_map_universe(isl_space_copy(p_space));
    for (size_t i = 0; i < top-1; ++i)
    {
      p_c = isl_constraint_alloc_equality(isl_local_space_copy(p_ls));
      p_c = isl_constraint_set_coefficient_si(p_c, isl_dim_out, i, 1);
      p_c = isl_constraint_set_coefficient_si(p_c, isl_dim_in, i, -1);
      p_tmp_map = isl_map_add_constraint(p_tmp_map, p_c);
    }

    p_c = isl_constraint_alloc_equality(isl_local_space_copy(p_ls));
    p_c = isl_constraint_set_coefficient_si(p_c, isl_dim_out, top-1, 1);
    p_c = isl_constraint_set_coefficient_si(p_c, isl_dim_in, top-1, -1);
    p_c = isl_constraint_set_constant_si(p_c, 1);
    p_tmp_map = isl_map_add_constraint(p_tmp_map, p_c);

    p_map = isl_map_union(p_map, p_tmp_map);
  }

  if (top < n_in_dims)
  {
    auto p_tmp_map = isl_map_lex_gt(isl_space_set_alloc(
      GetIslCtx().get(),
      isl_map_dim(p_map, isl_dim_param),
      n_in_dims - top
    ));
    p_tmp_map = insert_equal_dims(p_tmp_map, 0, 0, top);

    p_map = isl_map_union(p_map, p_tmp_map);
  }

  isl_local_space_free(p_ls);

  return isl::manage(p_map);
}

};  // namespace isl
