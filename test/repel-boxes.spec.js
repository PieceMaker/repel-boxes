const chai = require('chai');
const expect = chai.expect;
const RepelBoxes = require('../dist/repelBoxes');
const repelBoxes = new RepelBoxes();

describe('_toGeoJSONPoint', () => {

    it('should convert the origin to GeoJSON point', () => {
        const originGeoJSON = {type: "Point", coordinates: [0, 0]};
        expect(repelBoxes._toGeoJSONPoint(0, 0)).eql(originGeoJSON);
    });

});
