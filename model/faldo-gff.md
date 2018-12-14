# implentation of Faldo ontology for GFF#
            
            *   < p1 > sio: has - value  "1" ^ ^ xsd: integer;
            *         sio: refers - to < exon1 >.
            *
            * < exon1 > rdf: type insdc: Exon;
            *           faldo: location < region1 >.
            *
            * < region1 > rdf: type  faldo: Region;
            *             faldo: begin < position1 >;
            *              faldo: end < position2 >.
            *
            * < position1 > rdf: type  faldo: ExactPosition, faldo: ForwardStrandPosition;
            *               faldo: position  2983;
            *               faldo: reference < chromosome >.
            *
            * < position2 > rdf: type  faldo: ExactPosition, faldo: ForwardStrandPosition;
            *               faldo: position   10815;
            *               faldo: reference < chromosome >.
