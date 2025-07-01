"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
exports.fi = fi;
exports.serialize = serialize;
exports.approximateRenderTime = approximateRenderTime;
function fi() {
    return undefined;
}
/**
 * Serialize data to a string representation
 */
function serialize(data) {
    if (typeof data === "object") {
        if (data === null)
            return "null";
        if (data instanceof Date)
            return `new Date("${data.toISOString()}")`;
        else if (Array.isArray(data)) {
            let out = "[";
            for (const elem of data) {
                const serialized = serialize(elem);
                if (serialized !== undefined)
                    out += `${serialized},`;
            }
            if (out !== "[")
                out = out.slice(0, -1);
            out += "]";
            return out;
        }
        else {
            let out = "{";
            for (const prop in data) {
                const serialized = serialize(data[prop]);
                if (serialized !== undefined)
                    out += `"${prop}": ${serialized},`;
            }
            if (out !== "{")
                out = out.slice(0, -1);
            out += "}";
            return out;
        }
    }
    else {
        if (data !== undefined) {
            if (typeof data === "string")
                return `"${data.replace('"', '\\"').replace("\\", "\\\\")}"`;
            else if (typeof data === "boolean" || typeof data === "number")
                return data.toString();
        }
    }
}
// Parameters obtained by linear regression (R^2=0.9999) on experimental data
function approximateRenderTime(sz) {
    return 13.34865 * (sz / 1_000_000) + 0.63184;
}
//# sourceMappingURL=utils.js.map