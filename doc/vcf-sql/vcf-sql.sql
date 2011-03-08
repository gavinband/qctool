



-- ---
-- Globals
-- ---

-- SET SQL_MODE="NO_AUTO_VALUE_ON_ZERO";
-- SET FOREIGN_KEY_CHECKS=0;

-- ---
-- Table 'Symbol'
-- Defines symbols used by the db
-- ---

DROP TABLE IF EXISTS `Symbol`;
		
CREATE TABLE `Symbol` (
  `id` INTEGER AUTO_INCREMENT DEFAULT NULL,
  `Parent` INTEGER DEFAULT NULL,
  `Name` MEDIUMTEXT NOT NULL DEFAULT 'NULL',
  `Description` MEDIUMTEXT DEFAULT NULL,
  PRIMARY KEY (`id`)
) COMMENT='Defines symbols used by the db';

-- ---
-- Table 'Meta'
-- Contains metadata for this call set
-- ---

DROP TABLE IF EXISTS `Meta`;
		
CREATE TABLE `Meta` (
  `id` INTEGER AUTO_INCREMENT DEFAULT NULL,
  `Parent.id` INTEGER DEFAULT NULL,
  `Type` INTEGER NOT NULL DEFAULT NULL,
  `Value` MEDIUMTEXT DEFAULT NULL,
  PRIMARY KEY (`id`)
) COMMENT='Contains metadata for this call set';

-- ---
-- Table 'Variant'
-- 
-- ---

DROP TABLE IF EXISTS `Variant`;
		
CREATE TABLE `Variant` (
  `id` INTEGER AUTO_INCREMENT DEFAULT NULL,
  `Reference` INTEGER DEFAULT NULL,
  `Chrom` MEDIUMTEXT DEFAULT NULL,
  `Position` INTEGER NOT NULL DEFAULT NULL,
  `Name` MEDIUMTEXT DEFAULT NULL,
  PRIMARY KEY (`id`)
);

-- ---
-- Table 'VariantAttribute'
-- 
-- ---

DROP TABLE IF EXISTS `VariantAttribute`;
		
CREATE TABLE `VariantAttribute` (
  `Variant.id` INTEGER AUTO_INCREMENT DEFAULT NULL,
  `Type` INTEGER NOT NULL DEFAULT NULL,
  `Value` MEDIUMTEXT DEFAULT NULL,
  PRIMARY KEY (`Variant.id`)
);

-- ---
-- Table 'VariantCall'
-- Holds information about variant calls
-- ---

DROP TABLE IF EXISTS `VariantCall`;
		
CREATE TABLE `VariantCall` (
  `Variant.id` INTEGER AUTO_INCREMENT DEFAULT NULL,
  `Type` INTEGER NOT NULL DEFAULT NULL,
  `StorageType` INTEGER NOT NULL DEFAULT NULL,
  `DataLength` INTEGER NOT NULL DEFAULT NULL,
  `Data` BLOB NOT NULL DEFAULT 'NULL',
  PRIMARY KEY (`Variant.id`)
) COMMENT='Holds information about variant calls';

-- ---
-- Foreign Keys 
-- ---

ALTER TABLE `Symbol` ADD FOREIGN KEY (Parent) REFERENCES `Symbol` (`id`);
ALTER TABLE `Meta` ADD FOREIGN KEY (Parent.id) REFERENCES `Meta` (`id`);
ALTER TABLE `Meta` ADD FOREIGN KEY (Type) REFERENCES `Symbol` (`id`);
ALTER TABLE `Variant` ADD FOREIGN KEY (Reference) REFERENCES `Symbol` (`id`);
ALTER TABLE `VariantAttribute` ADD FOREIGN KEY (Variant.id) REFERENCES `Variant` (`id`);
ALTER TABLE `VariantAttribute` ADD FOREIGN KEY (Type) REFERENCES `Symbol` (`id`);
ALTER TABLE `VariantCall` ADD FOREIGN KEY (Variant.id) REFERENCES `Variant` (`id`);
ALTER TABLE `VariantCall` ADD FOREIGN KEY (Type) REFERENCES `Symbol` (`id`);
ALTER TABLE `VariantCall` ADD FOREIGN KEY (StorageType) REFERENCES `Symbol` (`id`);

-- ---
-- Table Properties
-- ---

-- ALTER TABLE `Symbol` ENGINE=InnoDB DEFAULT CHARSET=utf8 COLLATE=utf8_bin;
-- ALTER TABLE `Meta` ENGINE=InnoDB DEFAULT CHARSET=utf8 COLLATE=utf8_bin;
-- ALTER TABLE `Variant` ENGINE=InnoDB DEFAULT CHARSET=utf8 COLLATE=utf8_bin;
-- ALTER TABLE `VariantAttribute` ENGINE=InnoDB DEFAULT CHARSET=utf8 COLLATE=utf8_bin;
-- ALTER TABLE `VariantCall` ENGINE=InnoDB DEFAULT CHARSET=utf8 COLLATE=utf8_bin;

-- ---
-- Test Data
-- ---

-- INSERT INTO `Symbol` (`id`,`Parent`,`Name`,`Description`) VALUES
-- ('','','','');
-- INSERT INTO `Meta` (`id`,`Parent.id`,`Type`,`Value`) VALUES
-- ('','','','');
-- INSERT INTO `Variant` (`id`,`Reference`,`Chrom`,`Position`,`Name`) VALUES
-- ('','','','','');
-- INSERT INTO `VariantAttribute` (`Variant.id`,`Type`,`Value`) VALUES
-- ('','','');
-- INSERT INTO `VariantCall` (`Variant.id`,`Type`,`StorageType`,`DataLength`,`Data`) VALUES
-- ('','','','','');


