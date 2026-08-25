c Generated matrix elements still call these two split-order entry points.

      integer function orders_to_amp_split_pos(ord)
      use process_dimensions, only: nsplitorders
      use split_orders, only: module_orders_to_pos =>
     &     orders_to_amp_split_pos
      implicit none
      integer ord(*)

      orders_to_amp_split_pos=
     &     module_orders_to_pos(ord(1:nsplitorders))
      return
      end


      logical function orders_equal(orders1,orders2)
      use process_dimensions, only: nsplitorders
      use split_orders, only: module_orders_equal => orders_equal
      implicit none
      integer orders1(*),orders2(*)

      orders_equal=module_orders_equal(orders1(1:nsplitorders),
     &     orders2(1:nsplitorders))
      return
      end
