
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include "genfile/Error.hpp"
#include "components/CallComparerComponent/CallComparerComponent.hpp"
#include "components/CallComparerComponent/CallComparerDBOutputter.hpp"

CallComparerDBOutputter::UniquePtr CallComparerDBOutputter::create( std::string const& filename, std::string const& analysis, Metadata const& metadata ) {
	return UniquePtr( new CallComparerDBOutputter( filename, analysis, metadata ) ) ;
}

CallComparerDBOutputter::SharedPtr CallComparerDBOutputter::create_shared( std::string const& filename, std::string const& analysis, Metadata const& metadata ) {
	return SharedPtr( new CallComparerDBOutputter( filename, analysis, metadata ) ) ;
}

CallComparerDBOutputter::CallComparerDBOutputter( std::string const& filename, std::string const& analysis, Metadata const& metadata ):
	qcdb::DBOutputter( filename, analysis, metadata ),
	m_max_transaction_count( 10000 ),
	m_callset_id( get_or_create_entity( "callset", "entity representing call sets" ))
{
	db::Connection::ScopedTransactionPtr transaction = connection().open_transaction( 30 ) ;
	connection().run_statement(
		"CREATE TABLE IF NOT EXISTS CallComparison ( "
		"variant_id INT, callset1_id INT, callset2_id INT, method_id INT, variable_id INT, value FLOAT, "
		"FOREIGN KEY( variant_id ) REFERENCES Variant( id ), "
		"FOREIGN KEY( callset1_id ) REFERENCES Entity( id ), "
		"FOREIGN KEY( callset2_id ) REFERENCES Entity( id ), "
		"FOREIGN KEY( method_id ) REFERENCES Entity( id ), "
		"FOREIGN KEY( variable_id ) REFERENCES Entity( id ))"
	) ;
	connection().run_statement(
		"CREATE INDEX IF NOT EXISTS CallComparisonIndex ON CallComparison( variant_id, method_id, variable_id )"
	) ;

	connection().run_statement(
		"CREATE VIEW IF NOT EXISTS CallComparisonView AS "
		"SELECT V.variant_id, V.chromosome, V.position, V.rsid, callset1_id, C1.name AS callset1, callset2_id, C2.name AS callset2, variable_id, Variable.name AS variable, value "
		"FROM CallComparison CC "
		"INNER JOIN Variant V "
		"      ON V.id = CC.variant_id "
		"INNER JOIN Entity C1 "
		"      ON C1.id = CC.callset1_id"
		"INNER JOIN Entity C2 "
		"      ON C2.id = CC.callset2_id"
		"INNER JOIN Entity Variable "
		"      ON Variable.id = CC.variable_id "
		";"
	) ;

	construct_statements() ;
}

CallComparerDBOutputter::~CallComparerDBOutputter() {
	write_data( m_data ) ;
}

void CallComparerDBOutputter::begin_comparisons( genfile::SNPIdentifyingData const& snp ) {
	m_snp = snp ;
}

void CallComparerDBOutputter::end_comparisons() {}

void CallComparerDBOutputter::set_result(
	std::string const& callset1,
	std::string const& callset2,
	std::string const& comparison_method,
	std::string const& comparison_variable,
	genfile::VariantEntry const& value
) {
	m_data.resize( m_data.size() + 1 ) ;
	m_data.back().get<0>() = m_snp ;
	m_data.back().get<1>() = callset1 ;
	m_data.back().get<2>() = callset2 ;
	m_data.back().get<3>() = comparison_method ;
	m_data.back().get<4>() = comparison_variable ;
	m_data.back().get<5>() = value ;

	if( m_data.size() == m_max_transaction_count ) {
		write_data( m_data ) ;
		m_data.clear() ;
	}
}

void CallComparerDBOutputter::set_result(
	std::string const& comparison_method,
	std::string const& comparison_variable,
	genfile::VariantEntry const& value
) {
	m_data.resize( m_data.size() + 1 ) ;
	m_data.back().get<0>() = m_snp ;
	m_data.back().get<1>() = "" ;
	m_data.back().get<2>() = "" ;
	m_data.back().get<3>() = comparison_method ;
	m_data.back().get<4>() = comparison_variable ;
	m_data.back().get<5>() = value ;

	if( m_data.size() == m_max_transaction_count ) {
		write_data( m_data ) ;
		m_data.clear() ;
	}
}

void CallComparerDBOutputter::construct_statements() {
	m_insert_comparison_statement = connection().get_statement(
		"INSERT INTO CallComparison ( variant_id, callset1_id, callset2_id, method_id, variable_id, value ) "
		"VALUES( ?1, ?2, ?3, ?4, ?5, ?6 )"
	) ;
}

void CallComparerDBOutputter::write_data( Data const& data ) {
	db::Connection::ScopedTransactionPtr transaction = connection().open_transaction( 240 ) ; // wait up to 4 minutes.

	if( !transaction.get() ) {
		throw genfile::OperationFailedError( "CallComparerComponent::write_data()", connection().get_spec(), "Opening transaction." ) ;
	}
	for( std::size_t i = 0; i < m_data.size(); ++i ) {
		store_comparison(
			data[i].get<0>(),
			data[i].get<1>(),
			data[i].get<2>(),
			data[i].get<3>(),
			data[i].get<4>(),
			data[i].get<5>()
		) ;
	}
}

void CallComparerDBOutputter::store_comparison(
	genfile::SNPIdentifyingData const& snp,
	std::string const& callset1,
	std::string const& callset2,
	std::string const& comparison_method,
	std::string const& comparison_variable,
	genfile::VariantEntry const& value
) {
	db::Connection::RowId const snp_id = get_or_create_variant( snp ) ;

	db::Connection::RowId method_id = get_or_create_entity( comparison_method, "Method for performing pairwise genotype call comparison" ) ;
	db::Connection::RowId callset1_id = get_or_create_entity( callset1, "A callset", m_callset_id ) ;
	db::Connection::RowId callset2_id = get_or_create_entity( callset2, "A callset", m_callset_id ) ;
	db::Connection::RowId variable_id = get_or_create_entity( comparison_variable, "\"" + comparison_variable + "\" value for pairwise genotype call comparison" ) ;

	if( callset1 != "" ) {
		m_insert_comparison_statement
			->bind( 1, snp_id )
			.bind( 2, callset1_id )
			.bind( 3, callset2_id )
			.bind( 4, method_id )
			.bind( 5, variable_id )
			.bind( 6, value  )
			.step()
		;
		m_insert_comparison_statement->reset() ;
	} else {
		insert_summary_data( snp_id, variable_id, value ) ;
	}
	
}
