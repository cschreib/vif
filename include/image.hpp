#ifndef IMAGE_HPP
#define IMAGE_HPP

#include "vec.hpp"
#include "math.hpp"

template<typename Type>
typename vec_t<2,Type>::effective_type enlarge(const vec_t<2,Type>& v, int_t pix,
    const typename vec_t<2,Type>::rtype& def = 0.0) {

    if (pix >= 0) {
        typename vec_t<2,Type>::effective_type r = dblarr(v.dims[0]+2*pix, v.dims[1]+2*pix) + def;

        for (int_t y : rgen(v.dims[0]))
        for (int_t x : rgen(v.dims[1])) {
            r(x+pix,y+pix) = v(x,y);
        }

        return r;
    } else {
        assert(uint_t(-2*pix) < v.dims[0] && uint_t(-2*pix) < v.dims[1]);

        typename vec_t<2,Type>::effective_type r = dblarr(v.dims[0]+2*pix, v.dims[1]+2*pix);

        for (int_t y : rgen(r.dims[0]))
        for (int_t x : rgen(r.dims[1])) {
            r(x,y) = v(x-pix,y-pix);
        }

        return r;
    }
}

// Get a sub region 'reg' inside an image. reg = {x0, y0, x1, y1}.
template<typename TypeV, typename TypeR = int_t>
void subregion(const vec_t<2,TypeV>& v, const vec_t<1,TypeR>& reg, vec1i& rr, vec1i& rs) {
    assert(n_elements(reg) == 4);
    int_t nvx = v.dims[0], nvy = v.dims[1];
    int_t nx = reg(2)-reg(0)+1, ny = reg(3)-reg(1)+1;

    vec1i vreg = reg;
    vec1i sreg = {0,0,nx-1,ny-1};

    if (reg(0) >= nvx || reg(0) > reg(2) || reg(1) >= nvy || reg(1) > reg(3)) {
        rr = rs = intarr(0);
        return;
    }

    if (reg(0) < 0) {
        vreg(0) = 0;
        sreg(0) += 0 - reg(0);
    }
    if (reg(2) >= nvx) {
        vreg(2) = nvx-1;
        sreg(2) -= reg(2) - (nvx-1);
    }

    if (reg(1) < 0) {
        vreg(1) = 0;
        sreg(1) += 0 - reg(1);
    }
    if (reg(3) >= nvy) {
        vreg(3) = nvy-1;
        sreg(3) -= reg(3) - (nvy-1);
    }

    vec1i vx = rgen(vreg(0), vreg(2));
    vec1i vy = rgen(vreg(1), vreg(3));
    vec1i sx = rgen(sreg(0), sreg(2));
    vec1i sy = rgen(sreg(1), sreg(3));

    rr = flatten(indgen(nvx,nvy)(vx,vy));
    rs = flatten(indgen(nx,ny)(sx,sy));
}

template<typename TypeV, typename TypeR>
typename vec_t<2,TypeV>::effective_type subregion(const vec_t<2,TypeV>& v,
    const vec_t<1,TypeR>& reg, const typename vec_t<2,TypeV>::rtype& def = 0.0) {

    vec1i rr, rs;
    subregion(v, reg, rr, rs);

    int_t nx = reg(2)-reg(0)+1, ny = reg(3)-reg(1)+1;
    vec_t<2,rtype_t<TypeV>> sub = replicate(rtype_t<TypeV>(def), nx, ny);

    sub[rs] = v[rr];

    return sub;
}

template<typename TypeV, typename TypeR>
typename vec_t<2,TypeV>::effective_type subregion(const vec_t<2,TypeV>& v,
    std::initializer_list<TypeR>&& reg, const typename vec_t<2,TypeV>::rtype& def = 0.0) {

    return subregion(v, vec_t<1,TypeR>(std::move(reg)), def);
}

template<typename TypeV, typename TypeD = double>
typename vec_t<2,TypeV>::effective_type translate(const vec_t<2,TypeV>& v, double dx, double dy,
    const typename vec_t<2,TypeV>::rtype& def = 0.0) {

    vec_t<2,rtype_t<TypeV>> trs = replicate(rtype_t<TypeV>(def), v.dims);
    for (uint_t x = 0; x < v.dims[0]; ++x)
    for (uint_t y = 0; y < v.dims[1]; ++y) {
        double tx = x - dx;
        double ty = y - dy;
        int_t rx = round(tx);
        int_t ry = round(ty);

        if (rx < 0 || rx > int_t(v.dims[0])-1 || ry < 0 || ry > int_t(v.dims[1])-1) continue;
        int_t nx, ny;
        if (tx - rx > 0) {
            nx = rx + 1;
            if (nx > int_t(v.dims[0]-1)) continue;
        } else {
            nx = rx;
            rx = nx - 1;
            if (rx < 0) continue;
        }
        if (ty - ry > 0) {
            ny = ry + 1;
            if (ny > int_t(v.dims[1]-1)) continue;
        } else {
            ny = ry;
            ry = ny - 1;
            if (ry < 0) continue;
        }

        trs(x,y) = v(rx,ry)*(nx - tx)*(ny - ty) +
                   v(nx,ry)*(tx - rx)*(ny - ty) +
                   v(rx,ny)*(nx - tx)*(ty - ry) +
                   v(nx,ny)*(tx - rx)*(ty - ry);
    }

    return trs;
}

vec2d circular_mask(vec1u dim, const vec1d& center, double radius) {
    if (n_elements(dim) == 0) {
        dim = {uint_t(radius), uint_t(radius)};
    } else if (n_elements(dim) == 1) {
        dim = {dim[0], dim[0]};
    }

    assert(n_elements(center) == 2);

    vec3d m = dblarr(dim(0), dim(1), 4);

    vec2d px = replicate(dindgen(dim(1)), dim(0));
    vec2d py = transpose(replicate(dindgen(dim(0)), dim(1)));

    m(_,_,0) = pow(px-0.5 - center(1),2) + pow(py-0.5 - center(0),2) <= radius*radius;
    m(_,_,1) = pow(px+0.5 - center(1),2) + pow(py-0.5 - center(0),2) <= radius*radius;
    m(_,_,2) = pow(px+0.5 - center(1),2) + pow(py+0.5 - center(0),2) <= radius*radius;
    m(_,_,3) = pow(px-0.5 - center(1),2) + pow(py+0.5 - center(0),2) <= radius*radius;

    return mean(m,2);
}

#endif
