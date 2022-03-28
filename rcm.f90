MODULE rcm

CONTAINS

subroutine i4_swap ( i, j )

!*****************************************************************************80
!
!! I4_SWAP swaps two I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) I, J.  On output, the values of I and
!    J have been interchanged.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k

  k = i
  i = j
  j = k

  return
end
subroutine i4vec_reverse ( n, a )

!*****************************************************************************80
!
!! I4VEC_REVERSE reverses the elements of an I4VEC.
!
!  Example:
!
!    Input:
!
!      N = 5,
!      A = ( 11, 12, 13, 14, 15 ).
!
!    Output:
!
!      A = ( 15, 14, 13, 12, 11 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input/output, integer ( kind = 4 ) A(N), the array to be reversed.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i

  do i = 1, n/2
    call i4_swap ( a(i), a(n+1-i) )
  end do

  return
end

subroutine degree ( root, adj_num, adj_row, adj, mask, deg, iccsze, ls, &
  node_num )

!*****************************************************************************80
!
!! DEGREE computes the degrees of the nodes in the connected component.
!
!  Discussion:
!
!    The connected component is specified by MASK and ROOT.
!    Nodes for which MASK is zero are ignored.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 January 2003
!
!  Author:
!
!    Original FORTRAN77 version by Alan George, Joseph Liu.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ROOT, the node that defines the connected 
!    component.
!
!    Input, integer ( kind = 4 ) ADJ_NUM, the number of adjacency entries.
!
!    Input, integer ( kind = 4 ) ADJ_ROW(NODE_NUM+1).  Information about 
!    row I is stored in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
!    Input, integer ( kind = 4 ) ADJ(ADJ_NUM), the adjacency structure.
!    For each row, it contains the column indices of the nonzero entries.
!
!    Input, integer ( kind = 4 ) MASK(NODE_NUM), is nonzero for those nodes 
!    which are to be considered.
!
!    Output, integer ( kind = 4 ) DEG(NODE_NUM), contains, for each  node in 
!    the connected component, its degree.
!
!    Output, integer ( kind = 4 ) ICCSIZE, the number of nodes in the 
!    connected component.
!
!    Output, integer ( kind = 4 ) LS(NODE_NUM), stores in entries 1 through 
!    ICCSIZE the nodes in the connected component, starting with ROOT, and 
!    proceeding by levels.
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
  implicit none

  integer ( kind = 4 ) adj_num
  integer ( kind = 4 ) node_num

  integer ( kind = 4 ) adj(adj_num)
  integer ( kind = 4 ) adj_row(node_num+1)
  integer ( kind = 4 ) deg(node_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iccsze
  integer ( kind = 4 ) ideg
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jstop
  integer ( kind = 4 ) jstrt
  integer ( kind = 4 ) lbegin
  integer ( kind = 4 ) ls(node_num)
  integer ( kind = 4 ) lvlend
  integer ( kind = 4 ) lvsize
  integer ( kind = 4 ) mask(node_num)
  integer ( kind = 4 ) nbr
  integer ( kind = 4 ) node
  integer ( kind = 4 ) root
!
!  The sign of ADJ_ROW(I) is used to indicate if node I has been considered.
!
  ls(1) = root
  adj_row(root) = -adj_row(root)
  lvlend = 0
  iccsze = 1
!
!  LBEGIN is the pointer to the beginning of the current level, and
!  LVLEND points to the end of this level.
!
  do

    lbegin = lvlend + 1
    lvlend = iccsze
!
!  Find the degrees of nodes in the current level,
!  and at the same time, generate the next level.
!
    do i = lbegin, lvlend

      node = ls(i)
      jstrt = -adj_row(node)
      jstop = abs ( adj_row(node+1) ) - 1
      ideg = 0

      do j = jstrt, jstop

        nbr = adj(j)

        if ( mask(nbr) /= 0 ) then

          ideg = ideg + 1

          if ( 0 <= adj_row(nbr) ) then
            adj_row(nbr) = -adj_row(nbr)
            iccsze = iccsze + 1
            ls(iccsze) = nbr
          end if

        end if

      end do

      deg(node) = ideg

    end do
!
!  Compute the current level width.
!
    lvsize = iccsze - lvlend
!
!  If the current level width is nonzero, generate another level.
!
    if ( lvsize == 0 ) then
      exit
    end if

  end do
!
!  Reset ADJ_ROW to its correct sign and return.
!
  do i = 1, iccsze
    node = ls(i)
    adj_row(node) = -adj_row(node)
  end do

  return
end

subroutine rcm ( root, adj_num, adj_row, adj, mask, perm, iccsze, node_num )

!*****************************************************************************80
!
!! RCM renumbers a connected component by the reverse Cuthill McKee algorithm.
!
!  Discussion:
!
!    The connected component is specified by a node ROOT and a mask.
!    The numbering starts at the root node.
!
!    An outline of the algorithm is as follows:
!
!    X(1) = ROOT.
!
!    for ( I = 1 to N-1 )
!      Find all unlabeled neighbors of X(I),
!      assign them the next available labels, in order of increasing degree.
!
!    When done, reverse the ordering.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 August 2013
!
!  Author:
!
!    Original FORTRAN77 version by Alan George, Joseph Liu.
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ROOT, the node that defines the connected
!    component.  It is used as the starting point for the RCM ordering.
!    1 <= ROOT <= NODE_NUM.
!
!    Input, integer ( kind = 4 ) ADJ_NUM, the number of adjacency entries.
!
!    Input, integer ( kind = 4 ) ADJ_ROW(NODE_NUM+1).  Information about 
!    row I is stored in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
!    Input, integer ( kind = 4 ) ADJ(ADJ_NUM), the adjacency structure.
!    For each row, it contains the column indices of the nonzero entries.
!
!    Input/output, integer ( kind = 4 ) MASK(NODE_NUM), a mask for the nodes.  
!    Only those nodes with nonzero input mask values are considered by the 
!    routine.  The nodes numbered by RCM will have their mask values 
!    set to zero.
!
!    Output, integer ( kind = 4 ) PERM(NODE_NUM), the RCM ordering.
!
!    Output, integer ( kind = 4 ) ICCSZE, the size of the connected component
!    that has been numbered.
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!    1 <= NODE_NUM.
!
!  Local Parameters:
!
!    Workspace, integer DEG(NODE_NUM), a temporary vector used to hold 
!    the degree of the nodes in the section graph specified by mask and root.
!
  implicit none

  integer ( kind = 4 ) adj_num
  integer ( kind = 4 ) node_num

  integer ( kind = 4 ) adj(adj_num)
  integer ( kind = 4 ) adj_row(node_num+1)
  integer ( kind = 4 ) deg(node_num)
  integer ( kind = 4 ) fnbr
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iccsze
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jstop
  integer ( kind = 4 ) jstrt
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lbegin
  integer ( kind = 4 ) lnbr
  integer ( kind = 4 ) lperm
  integer ( kind = 4 ) lvlend
  integer ( kind = 4 ) mask(node_num)
  integer ( kind = 4 ) nbr
  integer ( kind = 4 ) node
  integer ( kind = 4 ) perm(node_num)
  integer ( kind = 4 ) root
!
!  Make sure NODE_NUM is legal.
!
  if ( node_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RCM - Fatal error!'
    write ( *, '(a,i4)' ) '  Illegal input value of NODE_NUM = ', node_num
    write ( *, '(a,i4)' ) '  Acceptable values must be positive.'
    stop 1
  end if
!
!  Make sure ROOT is legal.
!
  if ( root < 1 .or. node_num < root ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RCM - Fatal error!'
    write ( *, '(a,i4)' ) '  Illegal input value of ROOT = ', root
    write ( *, '(a,i4)' ) '  Acceptable values are between 1 and ', node_num
    stop 1
  end if
!
!  Find the degrees of the nodes in the component specified by MASK and ROOT.
!
  call degree ( root, adj_num, adj_row, adj, mask, deg, iccsze, perm, node_num )

  mask(root) = 0

  if ( iccsze < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RCM - Fatal error!'
    write ( *, '(a,i4)' ) '  Inexplicable component size ICCSZE = ', iccsze
    stop 1
  end if
!
!  If the connected component is a singleton, there is no reordering to do.
!
  if ( iccsze == 1 ) then
    return
  end if
!
!  Carry out the reordering procedure.
!
!  LBEGIN and LVLEND point to the beginning and
!  the end of the current level respectively.
!
  lvlend = 0
  lnbr = 1

  do while ( lvlend < lnbr )

    lbegin = lvlend + 1
    lvlend = lnbr

    do i = lbegin, lvlend
!
!  For each node in the current level...
!
      node = perm(i)
      jstrt = adj_row(node)
      jstop = adj_row(node+1) - 1
!
!  Find the unnumbered neighbors of NODE.
!
!  FNBR and LNBR point to the first and last neighbors
!  of the current node in PERM.
!
      fnbr = lnbr + 1

      do j = jstrt, jstop

        nbr = adj(j)

        if ( mask(nbr) /= 0 ) then
          lnbr = lnbr + 1
          mask(nbr) = 0
          perm(lnbr) = nbr
        end if

      end do
!
!  If no neighbors, skip to next node in this level.
!
      if ( lnbr <= fnbr ) then
        cycle
      end if
!
!  Sort the neighbors of NODE in increasing order by degree.
!  Linear insertion is used.
!
      k = fnbr

      do while ( k < lnbr )

        l = k
        k = k + 1
        nbr = perm(k)

        do while ( fnbr < l )

          lperm = perm(l)

          if ( deg(lperm) <= deg(nbr) ) then
            exit
          end if

          perm(l+1) = lperm
          l = l - 1

        end do

        perm(l+1) = nbr

      end do

    end do

  end do
!
!  We now have the Cuthill-McKee ordering.  
!  Reverse it to get the Reverse Cuthill-McKee ordering.
!
  call i4vec_reverse ( iccsze, perm )

  return
end

END MODULE rcm
