namespace genfile {
	namespace impl {
		template< typename Num, unsigned int DPs >
		struct fixed_precision_policy: public boost::spirit::karma::real_policies< Num >
		{
			static bool trailing_zeros( Num n ) { return false ; }
			static int floatfield( Num n ) { return boost::spirit::karma::real_policies< Num >::fmtflags::fixed ; }
			static unsigned precision( Num n ) { return DPs ; }
		} ;

		typedef boost::spirit::karma::real_generator<double, fixed_precision_policy<double,5> > FormatFloat5dp ;
		typedef boost::spirit::karma::real_generator<double, fixed_precision_policy<double,6> > FormatFloat6dp ;
		typedef boost::spirit::karma::real_generator<double, fixed_precision_policy<double,10> > FormatFloat10dp ;
	}
}
