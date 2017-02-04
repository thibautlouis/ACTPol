module fortran
contains

subroutine computeBinMap(iy,ix,bin_ly,bin_lx,pMShift_11,pMShift_12,pMShift_22,cos_array,sin_array,bMap,bMap0,bMap_cos,bMap_sin,bMap_cos2,bMap_sin2,bMap_cossin,phlx,phly,type,trimAtL)
  implicit none
  integer,intent(in) :: iy(:),ix(:),type,trimAtL
  real(8),intent(in) :: bin_ly(:),bin_lx(:),phlx(:),phly(:),pMShift_11(:,:),pMShift_12(:,:),pMShift_22(:,:),cos_array(:,:),sin_array(:,:),bMap(:,:)
  real(8),intent(inout) :: bMap0(:,:),bMap_cos(:,:),bMap_sin(:,:),bMap_cos2(:,:),bMap_sin2(:,:),bMap_cossin(:,:)
  !Work variables
  integer :: i,px,py,idx,idy
  real(8) :: Lx,Ly,ang2,c,s,w
  integer :: pxs(size(phlx)), nx



  do i = 1,size(iy) 
     Ly=bin_ly(iy(i)+1)
     Lx=bin_lx(ix(i)+1)
     ang2=atan2(Ly,Lx)
     w=bMap(ix(i)+1,iy(i)+1)

     nx=0
     do px =1, size(phlx)
        if (phlx(px)<trimAtL+Lx.and. phlx(px)>-trimAtL+Lx) then
           nx=nx+1
           pxs(nx) = px
        end if
     end do

     idy=0
     do py = 1, size(phly)
        if (phly(py)<trimAtL+Ly .and. phly(py)>-trimAtL+Ly) then
           idy=idy+1
           do idx = 1, nx
              px = pxs(idx)
              c = cos_array(idx,idy)*cos(-2*ang2)-sin_array(idx,idy)*sin(-2*ang2)
              s = sin_array(idx,idy)*cos(-2*ang2)+cos_array(idx,idy)*sin(-2*ang2)

              bMap0(idx,idy)=bMap0(idx,idy)+w*pMShift_11(px,py)
              bMap_cos(idx,idy)=bMap_cos(idx,idy)+w*pMShift_12(px,py)*c
              bMap_sin(idx,idy)=bMap_sin(idx,idy)+w*pMShift_12(px,py)*s
              bMap_cos2(idx,idy)=bMap_cos2(idx,idy)+w*pMShift_22(px,py)*c**2
              bMap_sin2(idx,idy)=bMap_sin2(idx,idy)+w*pMShift_22(px,py)*s**2
              bMap_cossin(idx,idy)=bMap_cossin(idx,idy)+w*pMShift_22(px,py)*s*c

              
           end do
        end if
     end do
   end do
end subroutine
end module
