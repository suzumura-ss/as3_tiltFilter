package info.smoche 
{
	import flash.display.BitmapData;
	import info.smoche.TiltFilterBase.TiltFilter_apply;
	import info.smoche.TiltFilterBase.TiltFilter_free;
	import info.smoche.TiltFilterBase.TiltFilter_init;
	/**
	 * ...
	 * @author Toshiyuki Suzumura
	 */
	public class TiltFilter 
	{
		public function TiltFilter()
		{
		}
		
		public static function tilt(yaw:Number, pitch:Number, roll:Number, source:BitmapData):BitmapData
		{
			var b:BitmapData = source.clone();
			var f:int = TiltFilter_init(yaw, pitch, roll);
			try {
				TiltFilter_apply(f, b, 0, uint( -1), null);
			} catch (e:Object) {
				TiltFilter_free(f);
				throw e;
			}
			TiltFilter_free(f);	
			return b;
		}
	}
}
